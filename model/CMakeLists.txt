list( APPEND oops_mom5cice5_src_files
c_mom5cice5_obs_data.F90
c_mom5cice5_setup.F90
c_mom5cice5_delete.F90
c_mom5cice5_prepare_integration.F90
c_mom5cice5_propagate.F90
ErrorCovariance.cc
ErrorCovariance.h
Fields.cc
Fields.h
Fortran.h
Geometry.cc
Geometry.h
Gom.cc
Gom.h
Increment.cc
Increment.h
LinearObsOp.h
Loc.h
Model.cc
Model.h
ModelBias.h
ModelBiasIncrement.h
ModelBiasCovariance.h
Observation.cc
Observation.h
ObsBias.cc
ObsBias.h
ObsBiasCovariance.cc
ObsBiasCovariance.h
ObsBiasIncrement.cc
ObsBiasIncrement.h
ObsFraction.cc
ObsFraction.h
ObsFractionTLAD.cc
ObsFractionTLAD.h
ObsHelp.cc
ObsHelp.h
ObsSpace.cc
ObsSpace.h
ObsVec.cc
ObsVec.h
State.cc
State.h
Traits.h
Variables.cc
Variables.h
generate_locations.f90
propagate.f90
mom5cice5_configs.F90
mom5cice5_constants.f90
mom5cice5_thermo.F90
mom5cice5_fields.F90
mom5cice5_fields_interface.F90
mom5cice5_fraction_mod.F90
mom5cice5_geom_mod.F90
mom5cice5_goms_mod.F90
mom5cice5_locs_mod.F90
mom5cice5_ncutils_mod.F90
mom5cice5_obs_data.F90
mom5cice5_obsoper_mod.F90
mom5cice5_obs_vectors.F90
mom5cice5_invent_state.f90
mom5cice5_vars_mod.F90
mom5cice5_balanceop.F90
)

if ( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( cfdealloc C_F_NO_DEALLOC )
endif()

ecbuild_add_library( TARGET   mom5cice5
                     DEFINITIONS ${cfdealloc}
                     SOURCES  ${oops_mom5cice5_src_files}
                     LIBS     oops ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

