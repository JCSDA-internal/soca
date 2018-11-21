/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "src/Traits.h"
#include "oops/runs/Moments.h"
#include "src/Run/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::Moments<soca::Traits> moments;  
  run.execute(moments);
  return 0;
}
