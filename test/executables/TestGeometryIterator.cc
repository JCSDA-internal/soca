/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "soca/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/GeometryIterator.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::GeometryIterator<soca::Traits> tests;
  return run.execute(tests);
}

