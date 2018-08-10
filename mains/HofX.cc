/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "src/Traits.h"
#include "oops/runs/HofX.h"
#include "src/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  oops::HofX<soca::Traits> hofx;
  run.execute(hofx);
  return 0;
}
