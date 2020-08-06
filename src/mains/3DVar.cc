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
#include "ufo/ObsTraits.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"
#include "oops/generic/instantiateModelFactory.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "saber/oops/instantiateLocalizationFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateModelFactory<soca::Traits>();
  soca::instantiateBalanceOpFactory();
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  saber::instantiateLocalizationFactory<soca::Traits>();
  saber::instantiateCovarFactory<soca::Traits>();
  oops::Variational<soca::Traits, ufo::ObsTraits> var;
  return run.execute(var);
}
