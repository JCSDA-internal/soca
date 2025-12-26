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

#include "netcdf.h"  // NOLINT
#include "oops/util/Logger.h"
#include "torch/torch.h"
#include "torch/csrc/distributed/c10d/ProcessGroup.hpp"
#include "torch/csrc/distributed/c10d/ProcessGroupMPI.hpp"

#include "BaseEmul.h"
#include "IceNet.h"

// -----------------------------------------------------------------------------

namespace soca {

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
    IceEmul(const eckit::Configuration & config,
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
      // TODO(G): Definitely move this elsewhere and leverage soca and/or atlas io
      // Read additional config
      std::string pole;
      getConfig().get("domain.pole", pole);
      bool cleanData;
      getConfig().get("domain.clean data", cleanData);

      // Read the patterns/targets
      // Open the file
      int ncid;
      nc_open(fileName.c_str(), NC_NOWRITE, &ncid);

      // Get dimensions
      int varid;
      std::string varName = "sst_h";
      nc_inq_varid(ncid, varName.c_str(), &varid);
      int ndims;
      nc_inq_varndims(ncid, varid, &ndims);
      int dimids[3];
      size_t dimLen[3];
      size_t nj, ni;
      nc_inq_vardimid(ncid, varid, dimids);
      nc_inq_dimlen(ncid, dimids[1], &nj);
      nc_inq_dimlen(ncid, dimids[2], &ni);

      // Allocate the CICE variable
      std::vector<float> lat(nj * ni);
      std::vector<float> lon(nj * ni);
      std::vector<float> aice(nj * ni);
      std::vector<float> tsfc(nj * ni);
      std::vector<float> sst(nj * ni);
      std::vector<float> sss(nj * ni);
      std::vector<float> sice(nj * ni);
      std::vector<float> hi(nj * ni);
      std::vector<float> hs(nj * ni);
      std::vector<float> mask(nj * ni);
      std::vector<float> tair(nj * ni);

      // Read the CICE variables
      varName = "ULAT"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, lat.data());
      varName = "ULON"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, lon.data());
      varName = "aice_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, aice.data());
      varName = "Tsfc_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, tsfc.data());
      varName = "sst_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, sst.data());
      varName = "sss_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, sss.data());
      varName = "sice_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, sice.data());
      varName = "hi_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, hi.data());
      varName = "hs_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, hs.data());
      varName = "umask"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, mask.data());
      varName = "Tair_h"; nc_inq_varid(ncid, varName.c_str(), &varid);
      nc_get_var_float(ncid, varid, tair.data());
      nc_close(ncid);

      // Calculate the number of patterns per pe
      int localBatchSize(0);
      int numPatterns(0);
      for (size_t i = getComm().rank(); i < lat.size(); i += getComm().size()) {
        if (selectData(mask[i], lat[i], aice[i], cleanData, pole)) {
          localBatchSize+=1;
        }
        if (localBatchSize >= n) { break; }
      }
      MPI_Allreduce(&localBatchSize, &numPatterns,
                    1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

      oops::Log::info() << "Number of patterns: " << numPatterns << std::endl;

      torch::Tensor patterns = torch::empty({numPatterns, getInputSize()}, torch::kFloat32);
      torch::Tensor targets = torch::empty({numPatterns}, torch::kFloat32);
      std::vector<float> lat_out;
      std::vector<float> lon_out;
      int cnt(0);
      for (size_t i = getComm().rank(); i < lat.size(); i += getComm().size()) {
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
      torch::Tensor mean = global_sum / (numPatterns * static_cast<float>(getComm().size()));
      torch::Tensor std = torch::sqrt(global_sq_sum /
                 (numPatterns * static_cast<float>(getComm().size())) - torch::pow(mean, 2));

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
      torch::Tensor input = torch::ones({getInputSize()});
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
        for (size_t i = 0; i < getInputSize(); ++i) {
          input[i] = inputs[j][i];
        }

        // Run the input through the FFNN
        torch::Tensor prediction = getModel()->forward(input);

        // Store results
        ice_original.push_back(targets[j].item<float>());
        ice_ffnn.push_back(prediction.item<float>());

        // Compute the Jacobian
        torch::Tensor doutdx = getModel()->jac(input);
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
      int ncid;
      nc_create(fileNameResults.c_str(), NC_CLOBBER, &ncid);
      int dimid, dimid2;
      nc_def_dim(ncid, "n", ice_original.size(), &dimid);
      nc_def_dim(ncid, "n_inputs", getInputSize(), &dimid2);

      int lon_varid, lat_varid, aice_varid, aice_ffnn_varid;
      int dcdt_varid, dcds_varid, dcdhs_varid, dcdhi_varid;
      int dcdsi_varid, dcdtsfc_varid, dcdtair_varid;
      nc_def_var(ncid, "lon", NC_FLOAT, 1, &dimid, &lon_varid);
      nc_def_var(ncid, "lat", NC_FLOAT, 1, &dimid, &lat_varid);
      nc_def_var(ncid, "aice", NC_FLOAT, 1, &dimid, &aice_varid);
      nc_def_var(ncid, "aice_ffnn", NC_FLOAT, 1, &dimid, &aice_ffnn_varid);
      nc_def_var(ncid, "dcdt", NC_FLOAT, 1, &dimid, &dcdt_varid);
      nc_def_var(ncid, "dcds", NC_FLOAT, 1, &dimid, &dcds_varid);
      nc_def_var(ncid, "dcdhs", NC_FLOAT, 1, &dimid, &dcdhs_varid);
      nc_def_var(ncid, "dcdhi", NC_FLOAT, 1, &dimid, &dcdhi_varid);
      nc_def_var(ncid, "dcdsi", NC_FLOAT, 1, &dimid, &dcdsi_varid);
      nc_def_var(ncid, "dcdtsfc", NC_FLOAT, 1, &dimid, &dcdtsfc_varid);
      nc_def_var(ncid, "dcdtair", NC_FLOAT, 1, &dimid, &dcdtair_varid);

      nc_enddef(ncid);

      nc_put_var_float(ncid, lon_varid, lon.data());
      nc_put_var_float(ncid, lat_varid, lat.data());
      nc_put_var_float(ncid, aice_varid, ice_original.data());
      nc_put_var_float(ncid, aice_ffnn_varid, ice_ffnn.data());
      nc_put_var_float(ncid, dcdt_varid, dcdt.data());
      nc_put_var_float(ncid, dcds_varid, dcds.data());
      nc_put_var_float(ncid, dcdhs_varid, dcdhs.data());
      nc_put_var_float(ncid, dcdhi_varid, dcdhi.data());
      nc_put_var_float(ncid, dcdsi_varid, dcdsi.data());
      nc_put_var_float(ncid, dcdtsfc_varid, dcdtsfc.data());
      nc_put_var_float(ncid, dcdtair_varid, dcdtair.data());

      nc_close(ncid);
    }
  };
}  // namespace soca
