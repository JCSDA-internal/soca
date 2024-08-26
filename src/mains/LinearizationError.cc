/*
 * (C) Copyright 2024- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "oops/runs/Run.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/LinearizationError.h"
#include "soca/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateLinearModelFactory<soca::Traits>();
  oops::instantiateModelFactory<soca::Traits>();
  oops::LinearizationError<soca::Traits> linearizationError;
  return run.execute(linearizationError);
}
