/*
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "mains/CheckpointModel.h"
#include "soca/Traits.h"
#include "soca/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  soca::CheckpointModel checkpoint;
  run.execute(checkpoint);
  return 0;
}
