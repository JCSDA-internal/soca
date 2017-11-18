/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/Utils/random_f.h"
#include <random>

// -----------------------------------------------------------------------------
namespace soca {
// -----------------------------------------------------------------------------

void random_f(const int & nn, double * xx) {
  static std::mt19937 generator(7);
  static std::normal_distribution<double> distribution(0.0, 1.0);

  for (int jj = 0; jj < nn; ++jj) xx[jj] = distribution(generator);
}

// -----------------------------------------------------------------------------

}  // namespace soca
