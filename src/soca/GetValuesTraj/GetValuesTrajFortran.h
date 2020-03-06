/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_GETVALUESTRAJ_GETVALUESTRAJFORTRAN_H_
#define SOCA_GETVALUESTRAJ_GETVALUESTRAJFORTRAN_H_

#include "soca/Fortran.h"

namespace soca {

  extern "C" {
    void soca_getvaltraj_setup_f90(const F90getvaltraj &);
    void soca_getvaltraj_delete_f90(const F90getvaltraj &);
  }
}  // namespace soca
#endif  // SOCA_GETVALUESTRAJ_GETVALUESTRAJFORTRAN_H_
