/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "oops/runs/ConvertIncrement.h"
#include "oops/runs/Run.h"
#include "soca/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::ConvertIncrement<soca::Traits> convertincrement;
  return run.execute(convertincrement);
}
