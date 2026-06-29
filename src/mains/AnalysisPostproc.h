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
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/Increment.h"
#include "oops/base/InflationBase.h"
#include "oops/base/instantiateInflationFactory.h"
#include "oops/base/StateSet.h"
#include "oops/base/State.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "soca/Geometry/FmsInput.h"
#include "soca/Traits.h"
#include "soca/Utils/incrqc/include/soca_incr_qc.h"

namespace soca {

// -----------------------------------------------------------------------------
/// \brief Application to recenter an ensemble of SOCA backgrounds around a
/// deterministic state with optional ensemble inflation. The application reads
/// an ensemble of states and a deterministic trajectory, and computes:
/// - ensemble perturbations: ensemble members minus ensemble mean
/// - optional scalar or field-based inflation of ensemble perturbations
/// - ocean recentering increments for MOM6 IAU: deterministic state minus ensemble
///   mean plus inflation effects. Postprocessing includes appending vertical
///   geometry needed by MOM6 IAU, optional zeroing of specified variables and
///   optional precision adjustment. If there's no inflation, all ensemble members
///   see the same recentering increment.
/// - ensemble of recentered states for CICE restarts: ensemble member + recentering
///   increment, postprocessed with Soca2Cice variable change.
/// - saves final increments and optionally ensemble mean/variance statistics
/// The states are processed in parallel, with each MPI task handling a subset of
/// ensemble members.
class AnalysisPostproc : public oops::Application {
  typedef oops::Geometry<soca::Traits>       Geometry_;
  typedef oops::Increment<soca::Traits>      Increment_;
  typedef oops::Increment4D<soca::Traits>    Increment4D_;
  typedef oops::IncrementSet<soca::Traits>   IncrementSet_;
  typedef oops::InflationBase<soca::Traits>  InflationBase_;
  typedef oops::InflationFactory<soca::Traits> InflationFactory_;
  typedef oops::StateSet<soca::Traits>       StateSet_;
  typedef oops::VariableChange<soca::Traits> VariableChange_;

 public:
  // -----------------------------------------------------------------------------

  explicit AnalysisPostproc(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
      oops::instantiateInflationFactory<soca::Traits>();
    }
  static const std::string classname() {return "soca::AnalysisPostproc";}

  // -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const {
    // SOCA-specific: copy/create files used by FMS once, to avoid race condition
    // in FMS MPP when creating parallel geometries. As hacky as it gets.
    const eckit::LocalConfiguration geomConfig(fullConfig, "geometry");
    soca::FmsInput fmsInput(this->getComm(), geomConfig);
    fmsInput.updateNameList();
    // Get the MPI partition
    const size_t nmembers = fullConfig.getInt("nens");
    const size_t nlocalmembers = fullConfig.getInt("nens per MPI task", 1);
    if (nmembers % nlocalmembers != 0) {
      oops::Log::error() << "Number of ensemble members must be divisible by the number of "
                        << "ensemble members per MPI task." << std::endl;
      throw eckit::UserError("Invalid number of ensemble members", Here());
    }
    const size_t nsubens = nmembers / nlocalmembers;
    const int ntasks = this->getComm().size();
    const int mytask = this->getComm().rank();
    const int tasks_per_subens = ntasks / nsubens;
    // first member of my subensemble
    int mymember = mytask / tasks_per_subens + 1;
    oops::Log::info() << "Running " << nmembers << " ensemble members handled by "
                << ntasks << " total MPI tasks and "
                << nlocalmembers << " ensemble members per MPI task, with "
                << tasks_per_subens << " MPI tasks per subensemble." << std::endl;
    if (ntasks % nsubens != 0) {
      oops::Log::error() << "Number of MPI tasks must be divisible by the number of "
                         << "subensembles." << std::endl;
      throw eckit::UserError("Invalid number of MPI tasks", Here());
    }

    // Create communicator across state, for each subensemble
    std::string commNameStr = "comm_state_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commState = this->getComm().split(mymember, commName);
    const int subrank = commState.rank();
    oops::Log::debug() << " commState name = " << commState.name()
                       << ", rank = " << commState.rank()
                       << ", size = " << commState.size() << std::endl;

    // Create communicator across ensemble members, for each subdomain
    std::string commEnsNameStr = "comm_ens_" + std::to_string(subrank);
    char const *commEnsName = commEnsNameStr.c_str();
    eckit::mpi::Comm & commEns = this->getComm().split(subrank, commEnsName);
    oops::Log::debug() << " commEns name = " << commEns.name()
                       << ", rank = " << commEns.rank()
                       << ", size = " << commEns.size() << std::endl;

    // Setup the  geometry of the ensemble members
    const Geometry_ geometry(geomConfig, commState);
    // Read all background states in parallel
    const eckit::LocalConfiguration statesConfig(fullConfig, "backgrounds");
    StateSet_ ens(geometry, statesConfig, oops::mpi::myself(), commEns);
    // Compute background ensemble mean as a StateSet (for computing the recentering increment
    // as a difference between two States)
    StateSet_ ensMean = ens.ens_mean();
    oops::Log::info() << "Background ensemble mean: " << ensMean << std::endl;
    oops::Log::test() << "Background ensemble mean: " << ensMean << std::endl;

    // Copy the ensemble into an ensemble of increments (it's still states),
    // used for computing ensemble statistics and background perturbations
    oops::Variables socaIncrVars(fullConfig, "increment variables");
    IncrementSet_ incs(geometry, socaIncrVars, ens);
    oops::Log::info() << "Input background states: " << incs << std::endl;

    // Compute background ensemble stats, print, save if needed.
    // bgEnsMeanInc is an IncrementSet (for computing the ensemble increments as difference
    // between two Increments)
    auto[bgEnsMeanInc, bgEnsVar] = incs.ens_stats();
    if ( fullConfig.has("ensemble mean output") && (commEns.rank() == 0) ) {
      const eckit::LocalConfiguration ensMeanOutputConfig(fullConfig, "ensemble mean output");
      ensMean.write(ensMeanOutputConfig);
    }
    if ( fullConfig.has("ensemble variance output") && (commEns.rank() == 0) ) {
      const eckit::LocalConfiguration ensVarianceOutputConfig(fullConfig,
                                                              "ensemble variance output");
      bgEnsVar.write(ensVarianceOutputConfig);
    }
    // Compute the standard deviation from the variance in place
    for (size_t jt = 0; jt < bgEnsVar.local_time_size(); ++jt) {
      bgEnsVar[jt].sqrt();
    }
    oops::Log::info() << " Background ensemble standard deviation: " << bgEnsVar << std::endl;
    oops::Log::test() << " Background ensemble standard deviation: " << bgEnsVar << std::endl;
    // Compute background ensemble perturbations
    incs -= bgEnsMeanInc;
    oops::Log::info() << " Background ensemble perturbations: " << incs << std::endl;
    oops::Log::test() << " Background ensemble perturbations: " << incs << std::endl;

    // Use ensemble DA increments if available
    if (fullConfig.has("analysis increments")) {
      // Read analysis increments
      const eckit::LocalConfiguration analysisIncrConfig(fullConfig, "analysis increments");
      IncrementSet_ analysisIncrs(geometry, socaIncrVars, incs.times(), analysisIncrConfig,
                                  oops::mpi::myself(), commEns);
      oops::Log::info() << "Analysis increments: " << analysisIncrs << std::endl;
      oops::Log::test() << "Analysis increments: " << analysisIncrs << std::endl;
      // Inflate analysis ensemble if needed
      if (fullConfig.has("ensemble inflation")) {
        const eckit::LocalConfiguration inflConfig(fullConfig, "ensemble inflation");
        std::unique_ptr<InflationBase_> inflation(InflationFactory_::create(inflConfig,
                                                  geometry, ens, socaIncrVars));
        inflation->doInflation(analysisIncrs);
        oops::Log::info() << "Analysis increments after inflation: "
                          << analysisIncrs << std::endl;
        oops::Log::test() << "Analysis increments after inflation: "
                          << analysisIncrs << std::endl;
      }
      // Update ensemble mean
      ensMean += analysisIncrs.ens_mean();
      oops::Log::test() << "Analysis ensemble mean: " << ensMean << std::endl;
      oops::Log::info() << "Analysis ensemble mean: " << ensMean << std::endl;
      // For statistics, add increments to background perturbations
      incs += analysisIncrs;
      // Compute analysis ensemble stats, print, save if needed.
      IncrementSet_ anEnsVar = incs.ens_var();
      if ( fullConfig.has("analysis ensemble mean output") && (commEns.rank() == 0) ) {
        const eckit::LocalConfiguration ensMeanOutputConfig(fullConfig,
          "analysis ensemble mean output");
        ensMean.write(ensMeanOutputConfig);
      }
      if ( fullConfig.has("analysis ensemble variance output") && (commEns.rank() == 0) ) {
        const eckit::LocalConfiguration ensVarianceOutputConfig(fullConfig,
          "analysis ensemble variance output");
        anEnsVar.write(ensVarianceOutputConfig);
      }
      // Compute the standard deviation from the variance in place
      for (size_t jt = 0; jt < anEnsVar.local_time_size(); ++jt) {
        anEnsVar[jt].sqrt();
      }
      oops::Log::info() << "Analysis ensemble standard deviation: " << anEnsVar << std::endl;
      oops::Log::test() << "Analysis ensemble standard deviation: " << anEnsVar << std::endl;
      // Set increments to be analysis increments for the postprocessing below
      incs = analysisIncrs;
    // No analysis ensemble increments provided: apply inflation to background
    // ensemble perturbations directly (Inflation classes operate on analysis
    // increments, which are not available in this branch)
    } else if (fullConfig.has("ensemble inflation")) {
      // Save the original increments for differencing later
      IncrementSet_ origincs(incs);
      // Multiplicative inflation of the ensemble perturbations
      const eckit::LocalConfiguration inflConfig(fullConfig, "ensemble inflation");
      if (inflConfig.has("value")) {
        const double inflation = inflConfig.getDouble("value");
        incs *= inflation;
        oops::Log::info() << "Ensemble perturbations after scalar inflation :"
                          << incs << std::endl;
      }
      if (inflConfig.has("field")) {
        Increment_ weight(geometry, socaIncrVars, incs[0].validTime());
        const eckit::LocalConfiguration weightConf(inflConfig, "field");
        weight.read(weightConf);
        for (size_t jj = 0; jj < incs.size(); ++jj) {
          incs[jj].schur_product_with(weight);
        }
        oops::Log::info() << "Ensemble perturbations after field inflation :"
                          << incs << std::endl;
      }
      // Increments that need to be added to the ensemble backgrounds
      incs -= origincs;
      oops::Log::info() << "Increments after inflation: " << incs << std::endl;
      oops::Log::test() << "Increments after inflation: " << incs << std::endl;
    } else {
      // if there's no inflation and no analysis increments, all the increments to
      // ensemble backgrounds are zero
      incs.zero();
    }

    // Read the state to recenter around if needed
    if (fullConfig.has("recentering state")) {
      const eckit::LocalConfiguration centerConfig(fullConfig, "recentering state");
      StateSet_ centerState(geometry, centerConfig);

      // Compute the recentering increment as the difference between
      // the ensemble mean and the deterministic
      Increment4D_ recenteringIncr(geometry, socaIncrVars,
                                   centerState.times(), centerState.commTime());
      recenteringIncr.diff(centerState, ensMean);
      oops::Log::info() << "Recentering increment: " << recenteringIncr << std::endl;
      oops::Log::test() << "Recentering increment: " << recenteringIncr << std::endl;
      incs += recenteringIncr;
      oops::Log::info() << "Increments after inflation and recentering: " << incs << std::endl;
      oops::Log::test() << "Increments after inflation and recentering: " << incs << std::endl;
    }
    if (fullConfig.has("increment postprocessing")) {
      const eckit::LocalConfiguration incPostprocConfig(fullConfig, "increment postprocessing");
      postprocessIncrements(ens, incs, incPostprocConfig);
      oops::Log::info() << "Increments after inflation and recentering and postprocessing: "
                        << incs << std::endl;
      oops::Log::test() << "Increments after inflation and recentering and postprocessing: "
                        << incs << std::endl;
    }

    // Save the increments used to initialize the ensemble forecast
    eckit::LocalConfiguration outputIncrConfig(fullConfig, "output increments");
    incs.write(outputIncrConfig);

    // Postprocess analysis (for CICE restarts) if needed
    if (fullConfig.has("analysis postprocessing")) {
      for (size_t iens = 0; iens < incs.local_ens_size(); ++iens) {
        const size_t ensMember = incs.local_ens()[iens];
        oops::Log::info() << "updating ice state " << ensMember << ":" << ens[iens] << std::endl;
        for (size_t itime = 0; itime < ens.local_time_size(); ++itime) {
          ens(itime, iens) += incs(itime, iens);
        }
        oops::Log::info() << "updated ice state " << ensMember << ":" << ens[iens] << std::endl;
        // set up variable change
        eckit::LocalConfiguration varchangeConfig(fullConfig,
          "analysis postprocessing.sea ice variable change");
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
  void postprocessIncrements(StateSet_ & ens,
                             IncrementSet_ & incs,
                             const eckit::Configuration & incPostprocConfig) const {
    // Add vertical geometry for MOM6 IAU
    if (incPostprocConfig.has("append vertical geometry")) {
      const eckit::LocalConfiguration vertGeomConfig(
        incPostprocConfig, "append vertical geometry");
      const eckit::LocalConfiguration fileConfig(vertGeomConfig, "vertical geometry");
      const std::string layerVarName = vertGeomConfig.getString("layers variable");
      oops::Variables layerVar({layerVarName});
      Increment_ vertGeom(incs.geometry(), layerVar, incs[0].validTime());
      vertGeom.read(fileConfig);
      oops::Variables vars = incs.variables();
      vars += layerVar;
      // Update the recentering increment with the vertical geometry
      for (size_t jj = 0; jj < incs.size(); ++jj) {
        // Note: these are soca::Increment specific methods
        incs[jj].increment().updateFields(vars);
        incs[jj].increment().updateFields(vertGeom.increment());
      }
    }
    // Set some variables to zero if needed
    if (incPostprocConfig.has("set increment variables to zero")) {
      oops::Variables socaZeroIncrVar(incPostprocConfig, "set increment variables to zero");
      if (!(socaZeroIncrVar <= incs.variables())) {
        // Add variables that have to be set to zero to the increment variables,
        // so that the increment can be updated
        oops::Variables newvars = incs.variables();
        newvars += socaZeroIncrVar;
        for (size_t jj = 0; jj < incs.size(); ++jj) {
          incs[jj].increment().updateFields(newvars);
        }
      }
      for (size_t jj = 0; jj < incs.size(); ++jj) {
        // Note: this is soca::Increment specific method
        incs[jj].increment().zero(socaZeroIncrVar);
      }
    }
    // Cut to a custom precision if needed
    if (incPostprocConfig.has("change precision")) {
      const eckit::LocalConfiguration precConfig(incPostprocConfig, "change precision");
      const oops::Variables precVars(precConfig, "variables");
      if (!(precVars <= incs.variables())) {
        oops::Log::error() << "Variables to change precision must be a subset of increment "
                           << "variables" << std::endl;
        throw eckit::UserError("Invalid variables to change precision", Here());
      }
      const double precision = precConfig.getDouble("precision");
      for (size_t jj = 0; jj < incs.size(); ++jj) {
        for (const auto & var : precVars) {
          auto field = incs[jj].fieldSet().fieldSet()[var.name()];
          auto view = atlas::array::make_view<double, 2>(field);
          for (int jnode = 0; jnode < view.shape(0); ++jnode) {
            for (int jlevel = 0; jlevel < view.shape(1); ++jlevel) {
              view(jnode, jlevel) = std::round(view(jnode, jlevel) / precision) * precision;
            }
          }
        }
        incs[jj].synchronizeFields();
      }
    }
    // Keep the increments within specified bounds if needed
    if (incPostprocConfig.has("bounds check")) {
      const eckit::LocalConfiguration boundsConfig(incPostprocConfig, "bounds check");

      // Apply QC to each ensemble member's increment
      for (size_t iens = 0; iens < incs.local_ens_size(); ++iens) {
        for (size_t itime = 0; itime < incs.local_time_size(); ++itime) {
          oops::Log::info() << "Applying increment QC to ensemble member "
                            << incs.local_ens()[iens] << " at time " << itime << std::endl;
          soca::incrqc::qcIncrement(ens(itime, iens).state(), incs(itime, iens).increment(),
                                     boundsConfig, incs.geometry().geometry());
        }
      }
    }
  }
  // -----------------------------------------------------------------------------
  std::string appname() const {
    return "soca::AnalysisPostproc";
  }
};

}  // namespace soca
