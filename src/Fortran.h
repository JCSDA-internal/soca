/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SOCA_MODEL_SOCAFORTRAN_H_
#define SOCA_MODEL_SOCAFORTRAN_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace soca {
  
  // Geometry key type
  typedef int F90geom;
  // Model key type
  typedef int F90model;
  // Locations key type
  typedef int F90locs;
  // Goms key type
  typedef int F90goms;
  // Fields key type
  typedef int F90flds;
  // Trajectory key type
  typedef int F90traj;
  // Trajectory for getvalues key type
  typedef int F90getvaltraj;  
  // Background error covariance key type
  typedef int F90bmat;
  // Observation vector key type
  typedef int F90ovec;
  // Obs operator key type
  typedef int F90hop;
  // Observation data base type
  typedef int F90odb;
  // Localization matrix
  typedef int F90lclz;

  /// Interface to Fortran SOCA model
  /*!
   * The core of the SOCA model is coded in Fortran.
   * Here we define the interfaces to the Fortran code.
   */

  extern "C" {
    // -----------------------------------------------------------------------------
    //  Geometry
    // -----------------------------------------------------------------------------
    void soca_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
    void soca_geo_clone_f90(const F90geom &, F90geom &);
    void soca_geo_info_f90(const F90geom &);
    void soca_geo_delete_f90(F90geom &);


    // -----------------------------------------------------------------------------
    //  Run
    // -----------------------------------------------------------------------------
    void soca_setup_f(const eckit::Configuration * const *);
    void soca_finalize_f();
  
    // -----------------------------------------------------------------------------
    //  Model
    // -----------------------------------------------------------------------------

    void soca_setup_f90(const eckit::Configuration * const *, const F90geom &, F90model &);
    void soca_delete_f90(F90model &);

    void soca_prepare_integration_f90(const F90model &, const F90flds &);
    //void soca_prepare_integration_tl_f90(const int & keyConf, const int & keyFlds);
    //void soca_prepare_integration_ad_f90(const int & keyConf, const int & keyFlds);

    void soca_propagate_f90(const F90model &, const F90flds &);
    //void soca_prop_traj_f90(const int & keyConf, const int & keyFlds, int & keyTraj);
    //void soca_propagate_tl_f90(const int & keyConf, const int & keyFlds,
    //                         const int & keyTraj);
    //void soca_propagate_ad_f90(const int & keyConf, const int & keyFlds,
    //                         const int & keyTraj);
    //
    //void soca_wipe_traj_f90(int & keyTraj);
    //void soca_traj_minmaxrms_f90(const int & keyTraj, double &);
  
    // -----------------------------------------------------------------------------
    //  Fields
    // -----------------------------------------------------------------------------
    void soca_field_create_f90(F90flds &, const F90geom &,
			       const eckit::Configuration * const *);
    void soca_field_delete_f90(F90flds &);
    void soca_field_copy_f90(const F90flds &, const F90flds &);
    void soca_field_zero_f90(const F90flds & );
    void soca_field_self_add_f90(const F90flds &, const F90flds &);
    void soca_field_self_sub_f90(const F90flds &, const F90flds &);
    void soca_field_self_mul_f90(const F90flds &, const double &);
    void soca_field_axpy_f90(const F90flds &, const double &, const F90flds &);
    void soca_field_dot_prod_f90(const F90flds &, const F90flds &, double &);
    void soca_field_self_schur_f90(const F90flds &, const F90flds &);
    void soca_field_random_f90(const F90flds &);
    void soca_field_dirac_f90(const F90flds &, const eckit::Configuration * const *);
  
    void soca_field_add_incr_f90(const F90flds &, const F90flds &);
    void soca_field_diff_incr_f90(const F90flds &, const F90flds &, const F90flds &);

    void soca_field_change_resol_f90(const F90flds &, const F90flds &);

    void soca_field_read_file_f90(const F90flds &, const eckit::Configuration * const *,
				  util::DateTime * const *);
    void soca_field_write_file_f90(const F90flds &, const eckit::Configuration * const *,
				   const util::DateTime * const *);
    void soca_field_interp_f90(const F90flds &, const F90locs &, const eckit::Configuration * const *, const F90goms &);
    void soca_field_interp_tl_f90(const F90flds &, const F90locs &, const eckit::Configuration * const *,const F90goms &);
    void soca_field_interp_tl_traj_f90(const F90flds &, const F90locs &, const eckit::Configuration * const *,const F90goms &, const F90getvaltraj &);    
    void soca_field_interp_ad_f90(const F90flds &, const F90locs &, const eckit::Configuration * const *,const F90goms &, const F90getvaltraj &);
  
    void soca_field_ug_coord_f90(const F90flds &, const int &);
    void soca_field_field_to_ug_f90(const F90flds &, const int &);
    void soca_field_field_from_ug_f90(const F90flds &, const int &);
  
    void soca_field_gpnorm_f90(const F90flds &, const int &, double &);
    void soca_field_sizes_f90(const F90flds &, int &, int &, int &, int &, int &, int &);
    void soca_field_rms_f90(const F90flds &, double &);

    // -----------------------------------------------------------------------------
    //  Trajectory (&more) for interpolation
    // -----------------------------------------------------------------------------
    void soca_getvaltraj_setup_f90(const F90getvaltraj &);
    void soca_getvaltraj_delete_f90(const F90getvaltraj &);

    // -----------------------------------------------------------------------------
    //  Background error
    // -----------------------------------------------------------------------------
    void soca_b_setup_f90(F90bmat &, const eckit::Configuration * const *, const F90geom &, const F90flds &);
    void soca_b_delete_f90(F90bmat &);
    void soca_b_linearize_f90(const F90flds &, const F90geom &);
    void soca_b_mult_f90(const F90bmat &, const F90flds &, const F90flds &, const F90flds &);
    void soca_b_invmult_f90(const F90bmat &, const F90flds &, const F90flds &);
    void soca_b_randomize_f90(const F90bmat &, const F90flds &);

    // -----------------------------------------------------------------------------
    //  Localization matrix
    // -----------------------------------------------------------------------------
    void soca_localization_setup_f90(F90lclz &, const eckit::Configuration * const *,
				     const F90geom &);
    void soca_localization_delete_f90(F90lclz &);
    void soca_localization_mult_f90(const F90lclz &, const F90flds &);
    
  }
}  // namespace soca
#endif  // SOCA_MODEL_SOCAFORTRAN_H_
