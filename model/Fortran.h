
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

/// Interface to Fortran MOM5CICE5 model
/*!
 * The core of the MOM5CICE5 model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void mom5cice5_geo_setup_f90(int & keyGeom, const eckit::Configuration * const *);
  void mom5cice5_geo_clone_f90(const int & keyGeom, int & keyGeom_other);
  void mom5cice5_geo_info_f90(const int & keyGeom, int &, int &, int &, int &, int &);
  void mom5cice5_geo_delete_f90(int & keyGeom);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
  void mom5cice5_setup_f90(const eckit::Configuration * const *, const int &, int & keyConf);
  void mom5cice5_delete_f90(int &);

  void mom5cice5_prepare_integration_f90(const int & keyConf, const int & keyFlds);
  //void mom5cice5_prepare_integration_tl_f90(const int & keyConf, const int & keyFlds);
  //void mom5cice5_prepare_integration_ad_f90(const int & keyConf, const int & keyFlds);

  void mom5cice5_propagate_f90(const int & keyConf, const int & keyFlds);
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
  void mom5cice5_field_create_f90(int & keyFlds, const int &, const int & keyVars);
  void mom5cice5_field_delete_f90(int & keyFlds);

  void mom5cice5_field_copy_f90(const int & keyFlds, const int & keyFldsOther);
  void mom5cice5_field_zero_f90(const int & keyFlds);
  void mom5cice5_field_dirac_f90(const int & keyFlds, const eckit::Configuration * const *);
  void mom5cice5_field_self_add_f90(const int & keyFlds, const int & keyFldsOther);
  void mom5cice5_field_self_sub_f90(const int & keyFlds, const int & keyFldsOther);
  void mom5cice5_field_self_mul_f90(const int & keyFlds, const double &);
  void mom5cice5_field_axpy_f90(const int & keyFlds, const double &, const int & keyFldsOther);
  void mom5cice5_field_dot_prod_f90(const int & keyFlds, const int & keyFldsOther, double &);
  void mom5cice5_field_self_schur_f90(const int & keyFlds, const int & keyFldsOther);
  void mom5cice5_field_random_f90(const int & keyFlds);

  void mom5cice5_field_add_incr_f90(const int & keyFlds, const int & keyFldsOther);
  void mom5cice5_field_diff_incr_f90(const int & keyFlds, const int & keyFldsOther,
                              const int & keyFldsOther2);

  void mom5cice5_field_change_resol_f90(const int & keyFlds, const int & keyFldsOther);

  void mom5cice5_field_read_file_f90(const int & keyFlds, const eckit::Configuration * const *,
                              util::DateTime * const *);
  void mom5cice5_field_write_file_f90(const int & keyFlds, const eckit::Configuration * const *,
                               const util::DateTime * const *);
  void mom5cice5_field_interp_f90(const int & keyFlds, const int &,
                           const int & keyGoms);
  void mom5cice5_field_interp_tl_f90(const int & keyFlds, const int &,
                              const int & keyGoms);
  void mom5cice5_field_interp_ad_f90(const int & keyFlds, const int &,
                              const int & keyGoms);  
  void mom5cice5_field_convert_to_f90(const int &, const int &);
  void mom5cice5_field_convert_from_f90(const int &, const int &);
  void mom5cice5_field_gpnorm_f90(const int & keyFlds, const int &, double &);
  void mom5cice5_field_sizes_f90(const int & keyFlds, int &, int &, int &, int &, int &, int &);
  void mom5cice5_field_rms_f90(const int & keyFlds, double &);

// -----------------------------------------------------------------------------
//  Variables
// -----------------------------------------------------------------------------
  void mom5cice5_var_create_f90(int & keyVars, const eckit::Configuration * const *);
  void mom5cice5_var_clone_f90(const int & keyVars, int & keyVars_other);
  void mom5cice5_var_info_f90(const int & keyVars, int &);
  void mom5cice5_var_delete_f90(int & keyVars);


// -----------------------------------------------------------------------------
//  Locations
// -----------------------------------------------------------------------------
  void mom5cice5_loc_delete_f90(int & keyLoc);
  void mom5cice5_loc_nobs_f90(const int & , int &);

// -----------------------------------------------------------------------------
//  Local Values (GOM)
// -----------------------------------------------------------------------------
  void mom5cice5_gom_create_f90(int & keyGoms);//, const int &);
  void mom5cice5_gom_delete_f90(int & keyGoms);
  void mom5cice5_gom_zero_f90(const int & keyGoms);
  void mom5cice5_gom_dotprod_f90(const int & keyGoms, const int & keyGomsOther, double &);
  void mom5cice5_gom_minmaxavg_f90(const int & keyGoms, int &, double &, double &, double &);

// -----------------------------------------------------------------------------
//  Fraction observations
// -----------------------------------------------------------------------------
  void mom5cice5_fraction_setup_f90(int & keyOper, const eckit::Configuration * const *);
  void mom5cice5_fraction_delete_f90(int & keyOper);
  void mom5cice5_fraction_equiv_f90(const int & keyGoms, const int & keyOvec, const double &);
  void mom5cice5_fraction_equiv_tl_f90(const int & keyGoms, const int & keyOvec, const double &);
  void mom5cice5_fraction_equiv_ad_f90(const int & keyGoms, const int & keyOvec, const double &);

// -----------------------------------------------------------------------------
//  Observation Vectors
// -----------------------------------------------------------------------------
  void mom5cice5_obsvec_setup_f90(int & keyOvec, const int &, const int &);
  void mom5cice5_obsvec_clone_f90(const int & keyOvec, int & keyOvecOther);
  void mom5cice5_obsvec_delete_f90(int & keyOvec);

  void mom5cice5_obsvec_assign_f90(const int & keyOvec, const int & keyOvecOther);
  void mom5cice5_obsvec_zero_f90(const int & keyOvec);
  void mom5cice5_obsvec_mul_scal_f90(const int & keyOvec, const double &);
  void mom5cice5_obsvec_add_f90(const int & keyOvec, const int & keyOvecOther);
  void mom5cice5_obsvec_sub_f90(const int & keyOvec, const int & keyOvecOther);
  void mom5cice5_obsvec_mul_f90(const int & keyOvec, const int & keyOvecOther);
  void mom5cice5_obsvec_div_f90(const int & keyOvec, const int & keyOvecOther);
  void mom5cice5_obsvec_axpy_f90(const int & keyOvec, const double &, const int & keyOvecOther);
  void mom5cice5_obsvec_invert_f90(const int & keyOvec);
  void mom5cice5_obsvec_random_f90(const int & keyOvec);
  void mom5cice5_obsvec_dotprod_f90(const int & keyOvec, const int & keyOvecOther, double &);
  void mom5cice5_obsvec_minmaxavg_f90(const int & keyOvec, double &, double &, double &);
  void mom5cice5_obsvec_nobs_f90(const int & keyOvec, int &);

// -----------------------------------------------------------------------------
//  Observation Handler
// -----------------------------------------------------------------------------
  void mom5cice5_obsdb_setup_f90(int & keyHelp, const eckit::Configuration * const *);
  void mom5cice5_obsdb_delete_f90(int & keyHelp);
  void mom5cice5_obsdb_get_f90(const int & keyHelp, const int &, const char *,
                        const int &, const char *, const int & keyOvec);
  void mom5cice5_obsdb_put_f90(const int & keyHelp, const int &, const char *,
                        const int &, const char *, const int & keyOvec);
  void mom5cice5_obsdb_locations_f90(const int & keyHelp, const int &, const char *,
                              const util::DateTime * const *, const util::DateTime * const *,
                              int & keyLoc);
  void mom5cice5_obsdb_getgom_f90(const int & keyHelp, const int &, const char *,
                           const int & keyVars,
                           const util::DateTime * const *, const util::DateTime * const *,
                           int & keyGoms);
  void mom5cice5_obsdb_generate_f90(const int & keyHelp, const int &, const char *,
                             const eckit::Configuration * const *, const util::DateTime * const *,
                             const util::Duration * const *, const int &, int &);
  void mom5cice5_obsdb_seterr_f90(const int & keyHelp, const int & keyOper, const double &);
  void mom5cice5_obsdb_nobs_f90(const int & keyHelp, const int &, const char *, int &);
  void mom5cice5_obsoper_inputs_f90(const int & keyOper, int & keyVars);
};
// -----------------------------------------------------------------------------


}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_MOM5CICE5FORTRAN_H_
