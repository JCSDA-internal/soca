/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef SOCA_MODEL_RANDOM_F_H_
#define SOCA_MODEL_RANDOM_F_H_

namespace soca {
extern "C" {
  void random_f(const int &, double *);
}
}  // namespace soca

#endif  // SOCA_MODEL_RANDOM_F_H_
