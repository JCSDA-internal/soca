/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "soca/Traits.h"
#include "soca/Run/Run.h"
#include "test/interface/ObsErrorCovariance.h"

int main(int argc,  char ** argv) {
  soca::Run run(argc, argv);
  test::ObsErrorCovariance<soca::Traits> tests;
  run.execute(tests);
  return 0;
}

