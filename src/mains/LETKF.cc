/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"
#include "oops/runs/LocalEnsembleDA.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsErrorFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/instantiateObsLocFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsLocFactory<soca::Traits>();
  ufo::instantiateObsErrorFactory<ufo::ObsTraits>();
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::LocalEnsembleDA<soca::Traits, ufo::ObsTraits> letkf;
  return run.execute(letkf);
}
