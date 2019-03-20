/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "src/Traits.h"
#include "oops/runs/CheckpointModel.h"
#include "src/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  oops::CheckpointModel<soca::Traits> checkpoint;
  run.execute(checkpoint);
  return 0;
}
