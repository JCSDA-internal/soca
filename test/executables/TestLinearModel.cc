/*
* (C) Copyright 2025 NOAA/NWS/NCEP/EMC
*
* This software is licensed under the terms of the Apache Licence Version 2.0
* which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include "soca/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/LinearModel.h"
#include "saber/oops/instantiateCovarFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  saber::instantiateCovarFactory<soca::Traits>();
  test::LinearModel<soca::Traits> tests;
  return run.execute(tests);
}
