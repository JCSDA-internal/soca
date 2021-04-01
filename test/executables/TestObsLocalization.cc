/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"

#include "ioda/instantiateObsLocFactory.h"
#include "oops/runs/Run.h"
#include "test/interface/ObsLocalization.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ioda::instantiateObsLocFactory<soca::Traits, ufo::ObsTraits>();
  test::ObsLocalization<soca::Traits, ufo::ObsTraits> tests;
  return run.execute(tests);
}
