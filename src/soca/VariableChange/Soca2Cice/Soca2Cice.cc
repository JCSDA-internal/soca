/*
 * (C) Copyright 2022-2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/VariableChange/Soca2Cice/Soca2Cice.h"
#include "soca/VariableChange/Soca2Cice/Soca2CiceFortran.h"

namespace soca {

// -----------------------------------------------------------------------------

static VariableChangeMaker<Soca2Cice>
                             makerVariableChangeSoca2Cice_("Soca2Cice");

// -----------------------------------------------------------------------------

Soca2Cice::Soca2Cice(const Geometry & geom,
                             const eckit::Configuration & conf)
  : geom_(geom)
{
  const eckit::Configuration * soca2ciceconfig = &conf;
  soca_soca2cice_setup_f90(keySoca2Cice_,
                           &soca2ciceconfig,
                           geom_.toFortran());
}

// -----------------------------------------------------------------------------

Soca2Cice::~Soca2Cice() {}

// -----------------------------------------------------------------------------

void Soca2Cice::changeVar(const State & xin, State & xout) const
{
  xout = xin;
  xin.syncFromFieldset();
  xout.syncFromFieldset();
  soca_soca2cice_changevar_f90(keySoca2Cice_, geom_.toFortran(),
                               xin.toFortran(), xout.toFortran());
  xout.syncToFieldset();
}

// -----------------------------------------------------------------------------

void Soca2Cice::changeVarInverse(const State &, State &) const {
}

// -----------------------------------------------------------------------------

}  // namespace soca
