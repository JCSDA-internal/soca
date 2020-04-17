/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "soca/Traits.h"
#include "soca/Transforms/instantiateBalanceOpFactory.h"
#include "saber/oops/EstimateParams.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  soca::instantiateBalanceOpFactory();
  oops::Run run(argc, argv);
  saber::EstimateParams<soca::Traits> dir;
  return run.execute(dir);
}
