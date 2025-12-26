/*
 * (C) Copyright 2024- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "oops/runs/Run.h"
#include "oops/runs/GenHybridLinearModelCoeffs.h"
#include "soca/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::GenHybridLinearModelCoeffs<soca::Traits> genHybridLinearModelCoeffs;
  return run.execute(genHybridLinearModelCoeffs);
}
