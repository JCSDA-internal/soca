list( APPEND oops_mom5cice5_src_files
#Geometry
./Geometry/Geometry.cc
./Geometry/Geometry.h
./Geometry/mom5cice5_geom_mod.F90

#Fields
./Fields/Fields.cc
./Fields/Fields.h
./Fields/mom5cice5_fields.F90
./Fields/mom5cice5_fields_interface.F90
./Fields/mom5cice5_invent_state.f90

#State
./State/State.cc
./State/State.h

#Variables
./Variables/Variables.cc
./Variables/Variables.h
./Variables/mom5cice5_vars_mod.F90

#Increment
./Increment/Increment.cc
./Increment/Increment.h

#Model
./Model/Model.cc
./Model/Model.h
./Model/c_mom5cice5_setup.F90
./Model/c_mom5cice5_delete.F90
./Model/c_mom5cice5_prepare_integration.F90
./Model/c_mom5cice5_propagate.F90
./Model/propagate.f90
./Model/mom5cice5_configs.F90

#ModelAtLocations
./ModelAtLocations/Gom.cc
./ModelAtLocations/Gom.h
./ModelAtLocations/mom5cice5_goms_mod.F90

#ObsSpace
./ObsSpace/ObsSpace.cc
./ObsSpace/ObsSpace.h
./ObsSpace/ObsHelp.cc
./ObsSpace/ObsHelp.h
./ObsSpace/mom5cice5_obs_data.F90
./ObsSpace/c_mom5cice5_obs_data.F90
./ObsSpace/generate_locations.f90

#Locations
./Locations/Loc.h
./Locations/mom5cice5_locs_mod.F90

#ObsOperator
./ObsOperator/Observation.cc
./ObsOperator/Observation.h
./ObsOperator/LinearObsOp.h
./ObsOperator/ObsFraction.cc
./ObsOperator/ObsFraction.h
./ObsOperator/ObsFractionTLAD.cc
./ObsOperator/ObsFractionTLAD.h
./ObsOperator/mom5cice5_fraction_mod.F90
./ObsOperator/mom5cice5_obsoper_mod.F90

#Covariance
./Covariance/ErrorCovariance.cc
./Covariance/ErrorCovariance.h

#ObsVector
./ObsVector/ObsVec.cc
./ObsVector/ObsVec.h
./ObsVector/mom5cice5_obs_vectors.F90

#Utils
./Utils/mom5cice5_constants.f90
./Utils/mom5cice5_thermo.F90
./Utils/mom5cice5_ncutils_mod.F90
./Utils/mom5cice5_balanceop.F90

Traits.h
Fortran.h

ModelBias.h
ModelBiasIncrement.h
ModelBiasCovariance.h
ObsBias.cc
ObsBias.h
ObsBiasCovariance.cc
ObsBiasCovariance.h
ObsBiasIncrement.cc
ObsBiasIncrement.h






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

