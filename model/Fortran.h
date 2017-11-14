
#ifndef MOM5CICE5_MODEL_MOM5CICE5FORTRAN_H_
#define MOM5CICE5_MODEL_MOM5CICE5FORTRAN_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace mom5cice5 {
  
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

/// Interface to Fortran MOM5CICE5 model
/*!
 * The core of the MOM5CICE5 model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void mom5cice5_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
  void mom5cice5_geo_clone_f90(const F90geom &, F90geom &);
  void mom5cice5_geo_info_f90(const F90geom &, int &, int &, int &, int &, int &);
  void mom5cice5_geo_delete_f90(F90geom &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------

  void mom5cice5_setup_f90(const eckit::Configuration * const *, const F90geom &, F90model &);
  void mom5cice5_delete_f90(F90model &);

  void mom5cice5_prepare_integration_f90(const F90model &, const F90flds &);
  //void mom5cice5_prepare_integration_tl_f90(const int & keyConf, const int & keyFlds);
  //void mom5cice5_prepare_integration_ad_f90(const int & keyConf, const int & keyFlds);

  void mom5cice5_propagate_f90(const F90model &, const F90flds &);
  //void mom5cice5_prop_traj_f90(const int & keyConf, const int & keyFlds, int & keyTraj);
  //void mom5cice5_propagate_tl_f90(const int & keyConf, const int & keyFlds,
  //                         const int & keyTraj);
  //void mom5cice5_propagate_ad_f90(const int & keyConf, const int & keyFlds,
  //                         const int & keyTraj);
  //
  //void mom5cice5_wipe_traj_f90(int & keyTraj);
  //void mom5cice5_traj_minmaxrms_f90(const int & keyTraj, double &);
  
// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void mom5cice5_field_create_f90(F90flds &, const F90geom &, const F90vars &);
  void mom5cice5_field_delete_f90(F90flds &);

  void mom5cice5_field_copy_f90(const F90flds &, const F90flds &);
  void mom5cice5_field_zero_f90(const F90flds & );
  void mom5cice5_field_dirac_f90(const F90flds &, const eckit::Configuration * const *);
  void mom5cice5_field_self_add_f90(const F90flds &, const F90flds &);
  void mom5cice5_field_self_sub_f90(const F90flds &, const F90flds &);
  void mom5cice5_field_self_mul_f90(const F90flds &, const double &);
  void mom5cice5_field_axpy_f90(const F90flds &, const double &, const F90flds &);
  void mom5cice5_field_dot_prod_f90(const F90flds &, const F90flds &, double &);
  void mom5cice5_field_self_schur_f90(const F90flds &, const F90flds &);
  void mom5cice5_field_random_f90(const F90flds &);

  void mom5cice5_field_add_incr_f90(const F90flds &, const F90flds &);
  void mom5cice5_field_diff_incr_f90(const F90flds &, const F90flds &, const F90flds &);

  void mom5cice5_field_change_resol_f90(const F90flds &, const F90flds &);

  void mom5cice5_field_read_file_f90(const F90flds &, const eckit::Configuration * const *,
                              util::DateTime * const *);
  void mom5cice5_field_write_file_f90(const F90flds &, const eckit::Configuration * const *,
                               const util::DateTime * const *);
  void mom5cice5_field_interp_f90(const F90flds &, const F90locs &, const F90goms &);
  void mom5cice5_field_interp_tl_f90(const F90flds &, const F90locs &, const F90goms &);
  void mom5cice5_field_interp_ad_f90(const F90flds &, const F90locs &, const F90goms &);
  
  void mom5cice5_field_convert_to_f90(const F90flds &, const int &);
  void mom5cice5_field_convert_from_f90(const F90flds &, const int &);
  
  void mom5cice5_field_gpnorm_f90(const F90flds &, const int &, double &);
  void mom5cice5_field_sizes_f90(const F90flds &, int &, int &, int &, int &, int &, int &);
  void mom5cice5_field_rms_f90(const F90flds &, double &);

// -----------------------------------------------------------------------------
//  Variables
// -----------------------------------------------------------------------------
  void mom5cice5_var_create_f90(F90vars &, const eckit::Configuration * const *);
  void mom5cice5_var_clone_f90(const F90vars &, F90vars &);
  void mom5cice5_var_info_f90(const F90vars &, int &);
  void mom5cice5_var_delete_f90(F90vars &);

// -----------------------------------------------------------------------------
//  Locations
// -----------------------------------------------------------------------------
  void mom5cice5_loc_delete_f90(F90locs &);
  void mom5cice5_loc_nobs_f90(const F90locs & , int &);

// -----------------------------------------------------------------------------
//  Local Values (GOM or GeoVaLs)
// -----------------------------------------------------------------------------
  void mom5cice5_gom_create_f90(F90goms &);//, const int &);
  void mom5cice5_gom_delete_f90(F90goms &);
  void mom5cice5_gom_zero_f90(const F90goms &);
  void mom5cice5_gom_random_f90(const F90goms &);
  void mom5cice5_gom_mult_f90(const F90goms &, const double &);  
  void mom5cice5_gom_dotprod_f90(const F90goms &, const F90goms &, double &);
  void mom5cice5_gom_minmaxavg_f90(const F90goms &, int &, double &, double &, double &);
  void mom5cice5_gom_read_file_f90(const F90goms &, const eckit::Configuration * const *);
  void mom5cice5_gom_write_file_f90(const F90goms &, const eckit::Configuration * const *);
  
// -----------------------------------------------------------------------------
//  Fraction observations
// -----------------------------------------------------------------------------
  void mom5cice5_fraction_setup_f90(F90hop &, const eckit::Configuration * const *);
  void mom5cice5_fraction_delete_f90(F90hop &);
  void mom5cice5_fraction_equiv_f90(const F90goms &, const F90ovec &, const double &);
  void mom5cice5_fraction_equiv_tl_f90(const F90goms &, const F90ovec &, const double &);
  void mom5cice5_fraction_equiv_ad_f90(const F90goms &, const F90ovec &, const double &);

// -----------------------------------------------------------------------------
//  Observation Vectors
// -----------------------------------------------------------------------------
  void mom5cice5_obsvec_setup_f90(F90hop &, const int &, const int &);
  void mom5cice5_obsvec_clone_f90(const F90ovec &, F90ovec &);
  void mom5cice5_obsvec_delete_f90(F90ovec &);

  void mom5cice5_obsvec_assign_f90(const F90ovec &, const F90ovec &);
  void mom5cice5_obsvec_zero_f90(const F90ovec &);
  void mom5cice5_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void mom5cice5_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void mom5cice5_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void mom5cice5_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void mom5cice5_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void mom5cice5_obsvec_axpy_f90(const F90ovec &, const double &, const F90ovec &);
  void mom5cice5_obsvec_invert_f90(const F90ovec &);
  void mom5cice5_obsvec_random_f90(const F90ovec &);
  void mom5cice5_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void mom5cice5_obsvec_minmaxavg_f90(const F90ovec &, double &, double &, double &);
  void mom5cice5_obsvec_nobs_f90(const F90ovec &, int &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void mom5cice5_obsdb_setup_f90(F90odb &, const eckit::Configuration * const *);
  void mom5cice5_obsdb_delete_f90(F90odb &);
  void mom5cice5_obsdb_get_f90(const F90odb &, const int &, const char *,
                               const int &, const char *, const F90ovec &);
  void mom5cice5_obsdb_put_f90(const F90odb &, const int &, const char *,
                               const int &, const char *, const F90ovec &);
  void mom5cice5_obsdb_locations_f90(const F90odb &, const int &, const char *,
                              const util::DateTime * const *, const util::DateTime * const *,
                              F90locs &);
  void mom5cice5_obsdb_getgom_f90(const F90odb &, const int &, const char *, const F90vars &,
                           const util::DateTime * const *, const util::DateTime * const *,
                           F90goms &);
  void mom5cice5_obsdb_generate_f90(const F90odb &, const int &, const char *,
                             const eckit::Configuration * const *, const util::DateTime * const *,
                             const util::Duration * const *, const int &, int &);
  void mom5cice5_obsdb_seterr_f90(const F90odb &, const F90hop &, const double &);
  void mom5cice5_obsdb_nobs_f90(const F90odb &, const int &, const char *, int &);
  void mom5cice5_obsoper_inputs_f90(const F90hop &, F90vars &);
};
// -----------------------------------------------------------------------------


}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_MOM5CICE5FORTRAN_H_
