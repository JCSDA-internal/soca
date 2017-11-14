/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MOM5CICE5_MODEL_RANDOM_F_H_
#define MOM5CICE5_MODEL_RANDOM_F_H_

namespace mom5cice5 {
extern "C" {
  void random_f(const int &, double *);
}
}  // namespace mom5cice5

#endif  // MOM5CICE5_MODEL_RANDOM_F_H_
