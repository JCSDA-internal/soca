/*
* (C) Copyright 2024 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <mpi.h>
#include <thread>  // NOLINT
#include <tuple>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/filesystem/PathName.h"

#include "netcdf"
#include "oops/util/Logger.h"
#include "torch/torch.h"
#include "torch/csrc/distributed/c10d/ProcessGroup.hpp"
#include "torch/csrc/distributed/c10d/ProcessGroupMPI.hpp"

#include "BaseEmul.h"
#include "IceNet.h"

// -----------------------------------------------------------------------------

namespace soca {

  /// Utilities:
  std::vector<float> readCice(const std::string fileName, const std::string varName) {
    // Open CICE restart file for reading
    netCDF::NcFile ncFile(fileName, netCDF::NcFile::read);

    // Get dimensions
    int dimLon = ncFile.getDim("ni").getSize();
    int dimLat = ncFile.getDim("nj").getSize();

    // Read the CICE variable
    std::vector<float> iceField(dimLat * dimLon);
    ncFile.getVar(varName).getVar(iceField.data());

    return iceField;
  }

  // Check if data is in the domain
  bool selectData(const float mask, const float lat, const float aice,
                const bool cleanData, const std::string& pole) {
    if (pole == "north") {
        if (cleanData) {
            return mask == 1 && lat > 40.0 && aice > 0.0 && aice <= 1.0;
        } else {
            return lat > 60.0;
        }
    } else if (pole == "south") {
        if (cleanData) {
            return mask == 1 && lat < -40.0 && aice > 0.0 && aice <= 1.0;
        } else {
            return lat < -60.0;
        }
    } else {
        oops::Log::error() << "Invalid pole value: " << pole << std::endl;
        return false;
    }
  }

  // IceEmul class derived from BaseEmul
  class IceEmul : public BaseEmul<IceNet> {
   public:
    // Constructor
    explicit IceEmul(const eckit::Configuration & config,
                     const eckit::mpi::Comm & comm) : BaseEmul<IceNet>(config, comm) {}

    // -----------------------------------------------------------------------------
    // Override prepData in IceEmul
    std::tuple<torch::Tensor,
               torch::Tensor,
               std::vector<float>,
               std::vector<float>,
               torch::Tensor,
               torch::Tensor>
    prepData(const std::string& fileName, bool geoloc = false, int n = 400000) override {
      // TODO(G): Move this elsewhere and leverage soca and/or atlas io
      // Read additional config
      std::string pole;
      config_.get("domain.pole", pole);
      bool cleanData;
      config_.get("domain.clean data", cleanData);

      // Read the patterns/targets
      std::vector<float> lat = readCice(fileName, "ULAT");
      std::vector<float> lon = readCice(fileName, "ULON");
      std::vector<float> aice = readCice(fileName, "aice_h");
      std::vector<float> tsfc = readCice(fileName, "Tsfc_h");
      std::vector<float> sst = readCice(fileName, "sst_h");
      std::vector<float> sss = readCice(fileName, "sss_h");
      std::vector<float> sice = readCice(fileName, "sice_h");
      std::vector<float> hi = readCice(fileName, "hi_h");
      std::vector<float> hs = readCice(fileName, "hs_h");
      std::vector<float> mask = readCice(fileName, "umask");
      std::vector<float> tair = readCice(fileName, "Tair_h");

      // Calculate the number of patterns per pe
      int localBatchSize(0);
      int numPatterns(0);
      for (size_t i = comm_.rank(); i < lat.size(); i += comm_.size()) {
        if (selectData(mask[i], lat[i], aice[i], cleanData, pole)) {
          localBatchSize+=1;
        }
        if (localBatchSize >= n) { break; }
      }
      MPI_Allreduce(&localBatchSize, &numPatterns,
                    1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

      oops::Log::info() << "Number of patterns: " << numPatterns << std::endl;

      torch::Tensor patterns = torch::empty({numPatterns, inputSize_}, torch::kFloat32);
      torch::Tensor targets = torch::empty({numPatterns}, torch::kFloat32);
      std::vector<float> lat_out;
      std::vector<float> lon_out;
      int cnt(0);
      for (size_t i = comm_.rank(); i < lat.size(); i += comm_.size()) {
        if (selectData(mask[i], lat[i], aice[i], cleanData, pole)) {
          patterns[cnt][0] = tair[i];
          patterns[cnt][1] = tsfc[i];
          patterns[cnt][2] = sst[i];
          patterns[cnt][3] = sss[i];
          patterns[cnt][4] = hs[i];
          patterns[cnt][5] = hi[i];
          patterns[cnt][6] = sice[i];

          targets[cnt] = aice[i];
          lat_out.push_back(lat[i]);
          lon_out.push_back(lon[i]);
          cnt+=1;
        }
        if (cnt >= numPatterns) { break; }
      }

      // Compute local sum and sum of squares for mean and std calculation
      torch::Tensor local_sum = torch::sum(patterns, /*dim=*/0);
      torch::Tensor local_sq_sum = torch::sum(torch::pow(patterns, 2), /*dim=*/0);

      // Initialize tensors to store the global sum and sum of squares
      torch::Tensor global_sum = torch::zeros_like(local_sum);
      torch::Tensor global_sq_sum = torch::zeros_like(local_sq_sum);

      // Use MPI_Allreduce to sum up all local sums and sq_sums across all processes
      MPI_Allreduce(local_sum.data_ptr(), global_sum.data_ptr(),
                    global_sum.numel(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(local_sq_sum.data_ptr(), global_sq_sum.data_ptr(),
                    global_sq_sum.numel(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

      // Calculate the global mean and std deviation
      torch::Tensor mean = global_sum / (numPatterns * static_cast<float>(comm_.size()));
      torch::Tensor std = torch::sqrt(global_sq_sum /
                 (numPatterns * static_cast<float>(comm_.size())) - torch::pow(mean, 2));

      return std::make_tuple(patterns, targets, lon_out, lat_out, mean, std);
    }

    // -----------------------------------------------------------------------------
    // Override predict in IceEmul
    void predict(const std::string& fileName, const std::string& fileNameResults,
                 const int n = -999) override {
      // Read the inputs/targets
      auto result = prepData(fileName, true);
      torch::Tensor inputs = std::get<0>(result);
      torch::Tensor targets = std::get<1>(result);
      std::vector lon = std::get<2>(result);
      std::vector lat = std::get<3>(result);

      // Loop through the patterns and predict
      torch::Tensor input = torch::ones({inputSize_});
      std::vector<float> ice_original;
      std::vector<float> ice_ffnn;
      // TODO(G): Store the jacobian in a 2D array
      std::vector<float> dcdt;
      std::vector<float> dcds;
      std::vector<float> dcdtsfc;
      std::vector<float> dcdtair;
      std::vector<float> dcdhi;
      std::vector<float> dcdhs;
      std::vector<float> dcdsi;
      for (size_t j = 0; j < targets.size(0); ++j) {
        for (size_t i = 0; i < inputSize_; ++i) {
          input[i] = inputs[j][i];
        }

        // Run the input through the FFNN
        torch::Tensor prediction = model_->forward(input);

        // Store results
        ice_original.push_back(targets[j].item<float>());
        ice_ffnn.push_back(prediction.item<float>());

        // Compute the Jacobian
        torch::Tensor doutdx = model_->jac(input);
        // Save the Jacobian elements into individual arrays
        // TODO(G): Store the jacobian in a 2D array
        dcdtair.push_back(doutdx[0].item<float>());
        dcdtsfc.push_back(doutdx[1].item<float>());
        dcdt.push_back(doutdx[2].item<float>());
        dcds.push_back(doutdx[3].item<float>());
        dcdhs.push_back(doutdx[4].item<float>());
        dcdhi.push_back(doutdx[5].item<float>());
        dcdsi.push_back(doutdx[6].item<float>());
      }

      // Save the prediction and Jacobian
      // TODO(G): Move into a separate function
      netCDF::NcFile ncFile(fileNameResults, netCDF::NcFile::replace);
      netCDF::NcDim dim = ncFile.addDim("n", ice_original.size());
      netCDF::NcDim dim2 = ncFile.addDim("n_inputs", inputSize_);

      ncFile.addVar("lon", netCDF::ncFloat, dim).putVar(lon.data());
      ncFile.addVar("lat", netCDF::ncFloat, dim).putVar(lat.data());
      ncFile.addVar("aice", netCDF::ncFloat, dim).putVar(ice_original.data());
      ncFile.addVar("aice_ffnn", netCDF::ncFloat, dim).putVar(ice_ffnn.data());
      ncFile.addVar("dcdt", netCDF::ncFloat, {dim}).putVar(dcdt.data());
      ncFile.addVar("dcds", netCDF::ncFloat, {dim}).putVar(dcds.data());
      ncFile.addVar("dcdhs", netCDF::ncFloat, {dim}).putVar(dcdhs.data());
      ncFile.addVar("dcdhi", netCDF::ncFloat, {dim}).putVar(dcdhi.data());
      ncFile.addVar("dcdsi", netCDF::ncFloat, {dim}).putVar(dcdsi.data());
      ncFile.addVar("dcdtsfc", netCDF::ncFloat, {dim}).putVar(dcdtsfc.data());
      ncFile.addVar("dcdtair", netCDF::ncFloat, {dim}).putVar(dcdtair.data());
    }
  };
}  // namespace soca
