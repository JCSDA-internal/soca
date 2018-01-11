
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
// Variables key type
typedef int F90vars;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Fields key type
typedef int F90flds;
// Trajectory key type
typedef int F90traj;
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
  void soca_field_create_f90(F90flds &, const F90geom &, const F90vars *);
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
  void soca_field_interp_f90(const F90flds &, const F90locs &, const F90goms &);
  void soca_field_interp_tl_f90(const F90flds &, const F90locs &, const F90goms &);
  void soca_field_interp_ad_f90(const F90flds &, const F90locs &, const F90goms &);
  
  void soca_field_convert_to_f90(const F90flds &, const int &);
  void soca_field_convert_from_f90(const F90flds &, const int &);
  
  void soca_field_gpnorm_f90(const F90flds &, const int &, double &);
  void soca_field_sizes_f90(const F90flds &, int &, int &, int &, int &, int &, int &);
  void soca_field_rms_f90(const F90flds &, double &);

// -----------------------------------------------------------------------------
//  Variables
// -----------------------------------------------------------------------------
//  void soca_var_create_f90(F90vars &, const eckit::Configuration * const *);
//  void soca_var_clone_f90(const F90vars &, F90vars &);
//  void soca_var_info_f90(const F90vars &, int &);
//  void soca_var_delete_f90(F90vars &);

// -----------------------------------------------------------------------------
//  Locations
// -----------------------------------------------------------------------------
  void soca_loc_delete_f90(F90locs &);
  void soca_loc_nobs_f90(const F90locs & , int &);

// -----------------------------------------------------------------------------
//  Local Values (GOM or GeoVaLs)
// -----------------------------------------------------------------------------
  void soca_gom_setup_f90(F90goms &, const F90locs &, const F90vars *);  
  void soca_gom_create_f90(F90goms &);//, const int &);
  void soca_gom_delete_f90(F90goms &);
  void soca_gom_zero_f90(const F90goms &);
  void soca_gom_random_f90(const F90goms &);
  void soca_gom_mult_f90(const F90goms &, const double &);  
  void soca_gom_dotprod_f90(const F90goms &, const F90goms &, double &);
  void soca_gom_minmaxavg_f90(const F90goms &, int &, double &, double &, double &);
  void soca_gom_read_file_f90(const F90goms &, const eckit::Configuration * const *);
  void soca_gom_write_file_f90(const F90goms &, const eckit::Configuration * const *);
  
// -----------------------------------------------------------------------------
//  Fraction observations
// -----------------------------------------------------------------------------
  void soca_fraction_setup_f90(F90hop &, const eckit::Configuration * const *);
  void soca_fraction_delete_f90(F90hop &);
  void soca_fraction_equiv_f90(const F90goms &, const F90ovec &, const double &);
  void soca_fraction_equiv_tl_f90(const F90goms &, const F90ovec &, const double &);
  void soca_fraction_equiv_ad_f90(const F90goms &, const F90ovec &, const double &);

// -----------------------------------------------------------------------------
//  Observation Vectors
// -----------------------------------------------------------------------------
  void soca_obsvec_setup_f90(F90hop &, const int &, const int &);
  void soca_obsvec_clone_f90(const F90ovec &, F90ovec &);
  void soca_obsvec_delete_f90(F90ovec &);

  void soca_obsvec_assign_f90(const F90ovec &, const F90ovec &);
  void soca_obsvec_zero_f90(const F90ovec &);
  void soca_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void soca_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void soca_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void soca_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void soca_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void soca_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void soca_obsvec_invert_f90(const F90ovec &);
  void soca_obsvec_random_f90(const F90ovec &);
  void soca_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void soca_obsvec_minmaxavg_f90(const F90ovec &, double &, double &, double &);
  void soca_obsvec_nobs_f90(const F90ovec &, int &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void soca_obsdb_setup_f90(F90odb &, const eckit::Configuration * const *);
  void soca_obsdb_delete_f90(F90odb &);
  void soca_obsdb_get_f90(const F90odb &, const int &, const char *,
                               const int &, const char *, const F90ovec &);
  void soca_obsdb_put_f90(const F90odb &, const int &, const char *,
                               const int &, const char *, const F90ovec &);
  void soca_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              const util::DateTime * const *, const util::DateTime * const *,
                              F90locs &);
  void soca_obsdb_getgom_f90(const F90odb &, const int &, const char *, const F90vars &,
                           const util::DateTime * const *, const util::DateTime * const *,
                           F90goms &);
  void soca_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration * const *, const util::DateTime * const *,
                             const util::Duration * const *, const int &, int &);
  void soca_obsdb_seterr_f90(const F90odb &, const F90hop &, const double &);
  void soca_obsdb_nobs_f90(const F90odb &, const int &, const char *, int &);
  void soca_obsoper_inputs_f90(const F90hop &, F90vars &);
};
// -----------------------------------------------------------------------------


}  // namespace soca
#endif  // SOCA_MODEL_SOCAFORTRAN_H_
