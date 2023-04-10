/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_FORTRAN_H_
#define SOCA_FORTRAN_H_

namespace soca {

  /// key type for soca_geom_mod::soca_geom
  typedef int F90geom;

  /// key type for soca_geom_iter_mod::soca_geom_iter
  typedef int F90iter;

  /// key type for soca_model_mod::soca_model
  typedef int F90model;

  /// key type for ufo_geovals_mod::ufo_geovals
  typedef int F90goms;

  /// key type for soca_fields_mod::soca_fields
  typedef int F90flds;

  /// key type for soca_getvalues_mod::soca_getvalues
  typedef int F90getval;

  /// key type for soca_covariance_mod::soca_cov
  typedef int F90bmat;

  /// key type for the various Fortran Transform classes.
  /**
   * \todo These are separate classes in separate registries, so they should
   * be separate key types
   */
  typedef int F90balopmat;
  typedef int F90varchange;
}  // namespace soca
#endif  // SOCA_FORTRAN_H_
