/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This  software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "oops/runs/ConvertToStructuredGrid.h"
#include "oops/runs/Run.h"
#include "soca/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::ConvertToStructuredGrid<soca::Traits> converttostructuredgrid;
  return run.execute(converttostructuredgrid);
}

