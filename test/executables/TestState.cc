/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Traits.h"
#include "src/Run/Run.h"
#include "test/interface/State.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  test::State<soca::Traits> tests;
  run.execute(tests);
  return 0;
}
