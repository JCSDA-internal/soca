/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"
#include "oops/runs/StaticBInit.h"
#include "soca/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  oops::StaticBInit<soca::Traits> bmat;
  run.execute(bmat);
  return 0;
}
