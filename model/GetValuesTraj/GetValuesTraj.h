/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_NOTHING_H_
#define SOCA_MODEL_NOTHING_H_

#include <ostream>

#include "oops/util/Printable.h"
#include "Fortran.h"

namespace soca {

// -----------------------------------------------------------------------------

class GetValuesTraj : public util::Printable {
 public:
  GetValuesTraj();
  ~GetValuesTraj();

  int & toFortran() {return keyGetValuesTraj_;}
  const int & toFortran() const {return keyGetValuesTraj_;}
  
 private:
  void print(std::ostream &) const {}
  F90getvaltraj keyGetValuesTraj_;
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_MODEL_NOTHING_H_
