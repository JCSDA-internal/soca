list( APPEND oops_soca_src_files
#Geometry
./Geometry/Geometry.cc
./Geometry/Geometry.h
./Geometry/soca_geom_mod.F90
./Geometry/soca_model_geom_type.F90

#Fields
./Fields/Fields.cc
./Fields/Fields.h
./Fields/soca_fields_mod.F90
./Fields/soca_fields.interface.F90
./Fields/soca_interpfields_mod.F90

#State
./State/State.cc
./State/State.h

#Variables
#./Variables/Variables.cc
#./Variables/Variables.h
#./Variables/soca_vars_mod.F90

#Increment
./Increment/Increment.cc
./Increment/Increment.h

#Model
./Model/Model.cc
./Model/Model.h
./Model/soca_model.interface.F90
./Model/soca_model_mod.F90

#Model
./LinearModel/TlmId.cc
./LinearModel/TlmId.h

#Covariance
./Covariance/ErrorCovariance.cc
./Covariance/ErrorCovariance.h
./Covariance/soca_covariance.interface.F90
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
./Utils/soca_seaice_balanceop.F90
./Utils/random_vectors_mod.F90
./Utils/random_f.cc
./Utils/random_f.h
./Utils/soca_mom6sis2.F90
./Utils/soca_interp.F90
./Utils/soca_Butils.F90
./Utils/soca_bumpinterp2d_mod.F90

#GetValuesTraj
./GetValuesTraj/GetValuesTraj.cc
./GetValuesTraj/GetValuesTraj.h
./GetValuesTraj/soca_getvaltraj_mod.F90

#Balance Operators
./Transforms/Ktc/Ktc.cc
./Transforms/Ktc/Ktc.h
./Transforms/Ktc/soca_ktc.interface.F90
./Transforms/Ktc/soca_ktc_mod.F90
./Transforms/Kst/Kst.cc
./Transforms/Kst/Kst.h
./Transforms/Kst/soca_kst.interface.F90
./Transforms/Kst/soca_kst_mod.F90
./Transforms/Ksshts/Ksshts.cc
./Transforms/Ksshts/Ksshts.h
./Transforms/Ksshts/soca_ksshts.interface.F90
./Transforms/Ksshts/soca_ksshts_mod.F90
./Transforms/BkgErr/BkgErr.cc
./Transforms/BkgErr/BkgErr.h
./Transforms/BkgErr/soca_bkgerr.interface.F90
./Transforms/BkgErr/soca_bkgerr_mod.F90
./Transforms/instantiateBalanceOpFactory.h
./Transforms/VertConv/VertConv.cc
./Transforms/VertConv/VertConv.h
./Transforms/VertConv/soca_vertconv.interface.F90
./Transforms/VertConv/soca_vertconv_mod.F90

Traits.h
Fortran.h

ModelBias.h
ModelBiasIncrement.h
ModelBiasCovariance.h
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

