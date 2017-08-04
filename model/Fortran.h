
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
  void mom5cice5_geo_info_f90(const int & keyGeom, int &, int &);
  void mom5cice5_geo_delete_f90(int & keyGeom);

// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void mom5cice5_field_create_f90(int & keyFlds, const int &, const int & keyVars);
  void mom5cice5_field_delete_f90(int & keyFlds);

  void mom5cice5_field_copy_f90(const int & keyFlds, const int & keyFldsOther);
  void mom5cice5_field_zero_f90(const int & keyFlds);
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

  void mom5cice5_field_gpnorm_f90(const int & keyFlds, const int &, double &);
  void mom5cice5_field_sizes_f90(const int & keyFlds, int &, int &, int &, int &);
  void mom5cice5_field_rms_f90(const int & keyFlds, double &);

// -----------------------------------------------------------------------------
//  Variables
// -----------------------------------------------------------------------------
  void mom5cice5_var_create_f90(int & keyVars, const eckit::Configuration * const *);
  void mom5cice5_var_clone_f90(const int & keyVars, int & keyVars_other);
  void mom5cice5_var_info_f90(const int & keyVars, int &, int &);
  void mom5cice5_var_delete_f90(int & keyVars);
}

// -----------------------------------------------------------------------------

}  // namespace mom5cice5
#endif  // MOM5CICE5_MODEL_MOM5CICE5FORTRAN_H_
