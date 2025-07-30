/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/Increment.h"
#include "oops/base/StateSet.h"
#include "oops/base/State.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "soca/Traits.h"

namespace soca {

class ParallelEnsRecenter : public oops::Application {
  typedef oops::Geometry<soca::Traits>       Geometry_;
  typedef oops::IncrementSet<soca::Traits>   IncrementSet_;
  typedef oops::State<soca::Traits>          State_;
  typedef oops::StateSet<soca::Traits>       StateSet_;
  typedef oops::VariableChange<soca::Traits> VariableChange_;

 public:
  // -----------------------------------------------------------------------------

  explicit ParallelEnsRecenter(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
  static const std::string classname() {return "soca::ParallelEnsRecenter";}

  // -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    const size_t nmembers = fullConfig.getInt("nens");

    // Get the MPI partition: communicators for ensemble members
    const int ntasks = this->getComm().size();
    const int mytask = this->getComm().rank();
    const int tasks_per_member = ntasks / nmembers;
    int mymember = mytask / tasks_per_member + 1;

    oops::Log::info() << "Running " << nmembers << " ensemble members handled by "
                << ntasks << " total MPI tasks and "
                << tasks_per_member << " MPI tasks per member." << std::endl;
    ASSERT(ntasks%nmembers == 0);

    // Create communicator across state, for each member
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = this->getComm().split(mymember, commName);
    const int subrank = commMember.rank();
    oops::Log::info() << " commMember name = " << commMember.name()
                       << ", rank = " << commMember.rank()
                       << ", size = " << commMember.size() << std::endl;

    // Create communicator across ensemble members, for each subdomain
    std::string commStateNameStr = "comm_state_" + std::to_string(subrank);
    char const *commStateName = commStateNameStr.c_str();
    eckit::mpi::Comm & commState = this->getComm().split(subrank, commStateName);
    oops::Log::info() << " commState name = " << commState.name()
                       << ", rank = " << commState.rank()
                       << ", size = " << commState.size() << std::endl;

    // Setup the native geometry of the ens. members,
    // currently assumed to be the same as the deterministic
    const eckit::LocalConfiguration geomConfig(fullConfig, "geometry");
    oops::Log::info() << "geometry: " << std::endl << geomConfig << std::endl;
    const Geometry_ geometryNative(geomConfig, commMember);

    // Setup the ensemble B (and output) soca geometry
    oops::Log::info() << "====================== ens B geometry" << std::endl;
    const std::string outputGeometryKey = fullConfig.has("output geometry")
                      ? "output geometry"  // keep things backward compatible for now
                      : "geometry";        // and default to the input geometry
    const eckit::LocalConfiguration geomOutConfig(fullConfig, outputGeometryKey);
    const Geometry_ geometry(geomOutConfig, commMember);

    // Read all states in parallel
    const eckit::LocalConfiguration statesConfig(fullConfig, "soca ensemble");
    std::unique_ptr<StateSet_> ensNative =
         std::make_unique<StateSet_>(geometryNative, statesConfig, oops::mpi::myself(), commState);
    StateSet_ ens(geometry, *ensNative);
    ensNative.reset();  // free memory used by ensNative
    // compute ensemble mean
    StateSet_ ensMean = ens.ens_mean();
    oops::Log::test() << "Ensemble mean: " << ensMean << std::endl;
    // Copy into ensemble of increments (even though they are stored as increments, it's still
    // states), and free memory used by ens
    oops::Variables socaIncrVar(fullConfig, "increment variables");
    IncrementSet_ incs(geometry, socaIncrVar, ens, true);
    oops::Log::info() << "Input states: " << incs << std::endl;
    // Compute ensemble stats, print, save if needed
    auto[ensMeanInc, ensVar] = incs.ens_stats();
    oops::Log::info() << " Ensemble mean: " << ensMean << std::endl;
    oops::Log::info() << " Ensemble mean from incs: " << ensMeanInc << std::endl;
    if ( fullConfig.has("ensemble mean output") ) {
      const eckit::LocalConfiguration ensMeanOutputConfig(fullConfig, "ensemble mean output");
      ensMean.write(ensMeanOutputConfig);
    }
    if ( fullConfig.has("ensemble variance output") ) {
      const eckit::LocalConfiguration ensVarianceOutputConfig(fullConfig,
                                                              "ensemble variance output");
      ensVar.write(ensVarianceOutputConfig);
    }
    // Compute the standard deviation from the variance in place
    for (size_t jt = 0; jt < ensVar.local_time_size(); ++jt) {
      ensVar[jt].sqrt();
    }
    oops::Log::info() << " Ensemble standard deviation: " << ensVar << std::endl;
    // Demean the increments
    incs -= ensMeanInc;
    oops::Log::info() << " Demeaned increments: " << incs << std::endl;
    oops::Log::test() << " Demeaned increments: " << incs << std::endl;

    // Initialize the trajectory for recentering, interpolate to the output geometry
    const eckit::LocalConfiguration trajConfig(fullConfig, "trajectory");
    StateSet_ determTrajNative(geometryNative, trajConfig);
    StateSet_ determTraj(geometry, determTrajNative);

    // Compute the recentering increment as the difference between
    // the ensemble mean and the deterministic
    IncrementSet_ recenteringIncr(geometry, socaIncrVar, determTraj.times(), determTraj.commTime());
    recenteringIncr.diff(determTraj, ensMean);
    // TODO(AS) postProcIncr.setToZero(recenteringIncr);
    oops::Log::info() << "Recentering increment: " << recenteringIncr << std::endl;

    bool seaiceRecenter = fullConfig.getBool("sea ice recenter", false);
    eckit::LocalConfiguration outputIncrConfig(fullConfig, "output increments");
    for (size_t iens = 0; iens < incs.local_ens_size(); ++iens) {
      const size_t ensMember = incs.local_ens()[iens];
      // make a copy of the recentering increment
      IncrementSet_ incr(recenteringIncr);

      // TODO(AS) Append the vertical geometry (for MOM6 IAU)
      //  soca::Increment mom6_incr = postProcIncr.appendLayer(incr);
      oops::Log::info() << "recentering incr " << iens << ":" << incr << std::endl;

      // TODO(AS) Set variables to zero if specified in the configuration
      //  postProcIncr.setToZero(incr);

      // TODO(AS) Optionally apply inflation
      /*
      if (fullConfig.has("ensemble inflation.value")) {
        const double inflation = fullConfig.getDouble("ensemble inflation.value");
        incr *= inflation;
        oops::Log::info() << "incr after scalar inflation " << iens << ":"
                          << incr << std::endl;
      }
      if (fullConfig.has("ensemble inflation.field")) {
        IncrementSet_ weight(geometry, incr.variables(), incr.validTimes(), incr.commTime());
        const eckit::LocalConfiguration weightConf(fullConfig, "ensemble inflation.field");
        weight.read(weightConf);
        incr.schur_product_with(weight);
        oops::Log::info() << "incr after field inflation " << iens << ":"
                          << incr << std::endl;
      }*/
      // Save the increments used to initialize the ensemble forecast
      outputIncrConfig.set("member", ensMember);
      incr.write(outputIncrConfig);

      // recenter ice if needed
      if (seaiceRecenter) {
        // read state
        eckit::LocalConfiguration ensmem_config(fullConfig, "sea ice analysis");
        std::string pattern;
        fullConfig.get("sea ice analysis.pattern", pattern);
        // replace templated string if necessary
        if (!pattern.empty()) {
          util::seekAndReplace(ensmem_config, pattern, ensMember, 0);
        }
        // assuming that this option is only used when recentering increment is on the
        // same geometry
        StateSet_ ens_an(geometry, ensmem_config);
        oops::Log::info() << "recentering ice state " << iens << ":" << ens_an << std::endl;
        ens_an += recenteringIncr;
        oops::Log::info() << "recentered ice state " << iens << ":" << ens_an << std::endl;
        // set up variable change
        eckit::LocalConfiguration varchange_config(fullConfig, "sea ice variable change");
        util::seekAndReplace(varchange_config, pattern, ensMember, 0);
        VariableChange_ varchange(varchange_config, geometry);
        oops::Variables varout(varchange_config, "output variables");
        // output happens inside soca2cice
        for (size_t jt = 0; jt < ens_an.local_time_size(); ++jt) {
          varchange.changeVar(ens_an[jt], varout);
        }
      }
    }
    return 0;
  }

 private:
  // -----------------------------------------------------------------------------
  std::string appname() const {
    return "soca::ParallelEnsRecenter";
  }
};

}  // namespace soca
