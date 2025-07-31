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

// -----------------------------------------------------------------------------
/// \brief Application to recenter an ensemble of SOCA backgrounds around a
/// deterministic state. The application reads an ensemble of states and a
/// deterministic trajectory, and computes:
/// - ocean recentering increment for MOM6 IAU: deterministic state minus ensemble mean,
///   with appended vertical geometry. All ensemble members see the same recentering
///   increment.
/// - ensemble of recentered states for CICE restarts: ensemble member + recentering
///   increment, postprocessed with Soca2Cice.
/// The states are processed in parallel, with each MPI task handling a single
/// ensemble member.
class ParallelEnsRecenter : public oops::Application {
  typedef oops::Geometry<soca::Traits>       Geometry_;
  typedef oops::Increment<soca::Traits>      Increment_;
  typedef oops::IncrementSet<soca::Traits>   IncrementSet_;
  typedef oops::StateSet<soca::Traits>       StateSet_;
  typedef oops::VariableChange<soca::Traits> VariableChange_;

 public:
  // -----------------------------------------------------------------------------

  explicit ParallelEnsRecenter(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
  static const std::string classname() {return "soca::ParallelEnsRecenter";}

  // -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    // Get the MPI partition: communicators for ensemble members
    const size_t nmembers = fullConfig.getInt("nens");
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

    // Setup the  geometry of the ensemble members
    const eckit::LocalConfiguration geomConfig(fullConfig, "geometry");
    const Geometry_ geometry(geomConfig, commMember);

    // Read all states in parallel
    const eckit::LocalConfiguration statesConfig(fullConfig, "soca ensemble");
    StateSet_ ens(geometry, statesConfig, oops::mpi::myself(), commState);
    // Compute ensemble mean as a StateSet (for computing the recentering increment as a difference
    // between two States)
    StateSet_ ensMean = ens.ens_mean();
    oops::Log::test() << "Ensemble mean: " << ensMean << std::endl;
    // Copy the ensemble into an ensemble of increments (it's still states)
    oops::Variables socaIncrVars(fullConfig, "increment variables");
    IncrementSet_ incs(geometry, socaIncrVars, ens);
    oops::Log::info() << "Input states: " << incs << std::endl;
    // Compute ensemble stats, print, save if needed.
    // ensMeanInc is an IncrementSet (for computing the ensemble increments as difference between
    // two Increments)
    auto[ensMeanInc, ensVar] = incs.ens_stats();
    oops::Log::info() << " Ensemble mean: " << ensMean << std::endl;
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
    // Compute ensemble perturbations
    incs -= ensMeanInc;
    oops::Log::info() << " Ensemble perturbations: " << incs << std::endl;
    oops::Log::test() << " Ensemble perturbations: " << incs << std::endl;

    // Initialize the trajectory for recentering, interpolate to the output geometry
    const eckit::LocalConfiguration trajConfig(fullConfig, "trajectory");
    StateSet_ determTraj(geometry, trajConfig);

    // Compute the recentering increment as the difference between
    // the ensemble mean and the deterministic
    IncrementSet_ recenteringIncr(geometry, socaIncrVars,
                                  determTraj.times(), determTraj.commTime());
    recenteringIncr.diff(determTraj, ensMean);
    oops::Log::info() << "Recentering increment: " << recenteringIncr << std::endl;
    oops::Log::test() << "Recentering increment: " << recenteringIncr << std::endl;

    // Add vertical geometry for MOM6 IAU
    const eckit::LocalConfiguration vertGeomConfig(fullConfig, "vertical geometry");
    const util::DateTime date(fullConfig.getString("date"));
    const std::string layerVarName = fullConfig.getString("layers variable");
    const oops::Variables layerVar({layerVarName});
    Increment_ vertGeom(geometry, layerVar, date);
    vertGeom.read(vertGeomConfig);
    socaIncrVars += layerVar;
    // Update the recentering increment with the vertical geometry
    for (size_t jt = 0; jt < recenteringIncr.local_time_size(); ++jt) {
      // Note: these are soca::Increment specific methods
      recenteringIncr[jt].increment().updateFields(socaIncrVars);
      recenteringIncr[jt].increment().updateFields(vertGeom.increment());
    }
    oops::Log::info() << "Recentering increment with vertical geometry: "
                      << recenteringIncr << std::endl;

    // TODO(AS) Add inflation

    bool seaiceRecenter = fullConfig.getBool("sea ice recenter", false);
    eckit::LocalConfiguration outputIncrConfig(fullConfig, "output increments");

    for (size_t iens = 0; iens < incs.local_ens_size(); ++iens) {
      const size_t ensMember = incs.local_ens()[iens];

      // Save the increments used to initialize the ensemble forecast
      // TODO(AS): fix to use different incs once inflation is implemented
      outputIncrConfig.set("member", ensMember);
      recenteringIncr.write(outputIncrConfig);

      // recenter ice if needed
      if (seaiceRecenter) {
        oops::Log::info() << "recentering ice state " << ensMember << ":" << ens[iens] << std::endl;
        for (size_t itime = 0; itime < ens.local_time_size(); ++itime) {
          ens(itime, iens) += recenteringIncr[itime];
        }
        oops::Log::info() << "recentered ice state " << ensMember << ":" << ens[iens] << std::endl;
        // set up variable change
        eckit::LocalConfiguration varchangeConfig(fullConfig, "sea ice variable change");
        std::string pattern = varchangeConfig.getString("pattern");
        util::seekAndReplace(varchangeConfig, pattern, ensMember+1, 0);
        VariableChange_ varchange(varchangeConfig, geometry);
        oops::Variables varout(varchangeConfig, "output variables");
        // output happens inside soca2cice
        for (size_t itime = 0; itime < ens.local_time_size(); ++itime) {
          varchange.changeVar(ens(itime, iens), varout);
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
