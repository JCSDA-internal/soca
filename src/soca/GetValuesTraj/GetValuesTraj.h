/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GETVALUESTRAJ_GETVALUESTRAJ_H_
#define SOCA_GETVALUESTRAJ_GETVALUESTRAJ_H_

#include <ostream>

#include "oops/util/Printable.h"

namespace soca {

// -----------------------------------------------------------------------------

class GetValuesTraj : public util::Printable {
 public:
  struct Ftn{};
  GetValuesTraj();
  ~GetValuesTraj();

  Ftn * & toFortran() {return ftn_;}
  Ftn * const & toFortran() const {return ftn_;}

 private:
  void print(std::ostream &) const {}
  Ftn * ftn_;
};

// -----------------------------------------------------------------------------

}  // namespace soca

#endif  // SOCA_GETVALUESTRAJ_GETVALUESTRAJ_H_
