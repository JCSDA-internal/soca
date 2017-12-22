list( APPEND oops_soca_src_files
#Geometry
./Geometry/Geometry.cc
./Geometry/Geometry.h
./Geometry/soca_geom_mod.F90
./Geometry/soca_model_geom_type.F90
#./Geometry/soca_geom_mom6.F90

#Fields
./Fields/Fields.cc
./Fields/Fields.h
./Fields/soca_fields.F90
./Fields/soca_fields_interface.F90
./Fields/soca_invent_state.F90

#State
./State/State.cc
./State/State.h

#Variables
./Variables/Variables.cc
./Variables/Variables.h
./Variables/soca_vars_mod.F90

#Increment
#./Increment/Increment.cc
#./Increment/Increment.h

#Model
./Model/Model.cc
./Model/Model.h
./Model/c_soca_setup.F90
./Model/c_soca_delete.F90
./Model/c_soca_prepare_integration.F90
./Model/c_soca_propagate.F90
./Model/propagate.F90
./Model/soca_configs.F90

#ModelAtLocations
./ModelAtLocations/Gom.cc
./ModelAtLocations/Gom.h
./ModelAtLocations/soca_goms_mod.F90

#ObsSpace
./ObsSpace/ObsSpace.cc
./ObsSpace/ObsSpace.h
./ObsSpace/ObsHelp.cc
./ObsSpace/ObsHelp.h
./ObsSpace/soca_obs_data.F90
./ObsSpace/c_soca_obs_data.F90
./ObsSpace/generate_locations.F90

#Locations
./Locations/Loc.h
./Locations/soca_locs_mod.F90

#ObsOperator
#./ObsOperator/Observation.cc
#./ObsOperator/Observation.h
#./ObsOperator/LinearObsOp.h
#./ObsOperator/ObsFraction.cc
#./ObsOperator/ObsFraction.h
#./ObsOperator/ObsFractionTLAD.cc
#./ObsOperator/ObsFractionTLAD.h
#./ObsOperator/soca_fraction_mod.F90
./ObsOperator/soca_obsoper_mod.F90

#Covariance
#./Covariance/ErrorCovariance.cc
#./Covariance/ErrorCovariance.h

#ObsVector
./ObsVector/ObsVec.cc
./ObsVector/ObsVec.h
./ObsVector/soca_obs_vectors.F90

#Run
./Run/Run.cc
./Run/Run.h
./Run/soca_setup.f90

#Utils
./Utils/soca_constants.F90
./Utils/soca_thermo.F90
./Utils/soca_ncutils_mod.F90
#./Utils/soca_balanceop.F90
./Utils/random_vectors_mod.F90
./Utils/random_f.cc
./Utils/random_f.h
#./Utils/linkedList_c.f
#./Utils/linkedList_i.f
./Utils/soca_mom6sis2.F90

Traits.h
Fortran.h

#ModelBias.h
#ModelBiasIncrement.h
#ModelBiasCovariance.h
#ObsBias.cc
#ObsBias.h
#ObsBiasCovariance.cc
#ObsBiasCovariance.h
#ObsBiasIncrement.cc
#ObsBiasIncrement.h
)

if ( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( cfdealloc C_F_NO_DEALLOC )
endif()

#set ( FMS_LIBRARIES $ENV{FMS_LIBRARIES} )
#set ( FMS_INCLUDE_DIRS $ENV{FMS_INCLUDE_DIRS} )

ecbuild_add_library( TARGET   soca
                     DEFINITIONS ${cfdealloc}
                     SOURCES  ${oops_soca_src_files}
                     LIBS     oops ioda ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${MOM6_LIBRARIES} ${FMS_LIBRARIES} ${MPI_LIBRARIES}
		     INCLUDES ${FMS_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

