/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"
#include "oops/runs/EnsHofX.h"
#include "soca/instantiateObsFilterFactory.h"
#include "soca/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  soca::instantiateObsFilterFactory();
  oops::EnsHofX<soca::Traits> hofx;
  run.execute(hofx);
  return 0;
}
