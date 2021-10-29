/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/EnsembleApplication.h"
#include "oops/runs/HofX4D.h"
#include "oops/runs/Run.h"
#include "soca/Traits.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::EnsembleApplication<oops::HofX4D <soca::Traits, ufo::ObsTraits> > hofx;
  return run.execute(hofx);
}
