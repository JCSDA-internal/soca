/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_RUN_RUN_H_
#define SRC_RUN_RUN_H_

#include "oops/runs/Run.h"

namespace soca {

  // -----------------------------------------------------------------------------

  class Run : public oops::Run {
   public:
    Run(int, char **);
    ~Run();
  };

  // -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SRC_RUN_RUN_H_
