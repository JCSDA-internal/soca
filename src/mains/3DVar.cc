/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "soca/Traits.h"
#include "soca/Transforms/instantiateBalanceOpFactory.h"
#include "soca/instantiateObsFilterFactory.h"
#include "oops/runs/Variational.h"
#include "soca/Run/Run.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  soca::instantiateBalanceOpFactory();
  soca::instantiateObsFilterFactory();
  oops::Variational<soca::Traits> var;
  run.execute(var);
  return 0;
}
