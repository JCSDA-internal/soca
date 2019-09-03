/*
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"
#include "soca/Run/Run.h"
#include "test/interface/Variables.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  test::Variables<soca::Traits> tests;
  run.execute(tests);
  return 0;
}

