list( APPEND oops_soca_src_files
#Geometry
./Geometry/Geometry.cc
./Geometry/Geometry.h
./Geometry/soca_geom_mod.F90
./Geometry/soca_model_geom_type.F90

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
./Increment/Increment.cc
./Increment/Increment.h

#Model
./Model/Model.cc
./Model/Model.h
./Model/c_soca_setup.F90
./Model/c_soca_delete.F90
./Model/c_soca_prepare_integration.F90
./Model/c_soca_propagate.F90
./Model/propagate.F90
./Model/soca_configs.F90

#Model
./LinearModel/TlmId.cc
./LinearModel/TlmId.h

#Covariance
./Covariance/ErrorCovariance.cc
./Covariance/ErrorCovariance.h
./Covariance/soca_bmatrix.interface.F90
./Covariance/soca_covariance_mod.F90

#LocalizationMatrix
./LocalizationMatrix/LocalizationMatrix.cc
./LocalizationMatrix/LocalizationMatrix.h
./LocalizationMatrix/instantiateLocalizationFactory.h
./LocalizationMatrix/soca_localization_matrix.F90

#Run
./Run/Run.cc
./Run/Run.h
./Run/soca_setup.f90

#Utils
./Utils/soca_utils.F90
./Utils/soca_constants.F90
./Utils/soca_thermo.F90
./Utils/soca_balanceop.F90
./Utils/random_vectors_mod.F90
./Utils/random_f.cc
./Utils/random_f.h
./Utils/soca_mom6sis2.F90
./Utils/soca_interp.F90
./Utils/soca_Butils.F90

#Run
./Nothing/Nothing.h

Traits.h
Fortran.h

ModelBias.h
ModelBiasIncrement.h
ModelBiasCovariance.h
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
                     LIBS     oops ioda ufo ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES} ${MOM6_LIBRARIES} ${FMS_LIBRARIES} ${MPI_LIBRARIES}
		     INCLUDES ${FMS_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS} ${MOM6_INCLUDE_DIRS}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

