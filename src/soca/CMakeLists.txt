# (C) Copyright 2017-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND soca_src_files
#Geometry
Geometry/Geometry.cc
Geometry/Geometry.h
Geometry/soca_geom_mod.F90
Geometry/soca_geom_interface_mod.F90

#Fields
Fields/Fields.cc
Fields/Fields.h
Fields/soca_fields_mod.F90
Fields/soca_fields_interface_mod.F90
Fields/soca_interpfields_mod.F90
Fields/soca_ocnsfc_mod.F90
Fields/soca_seaice_mod.F90
Fields/soca_fieldsutils_mod.F90

#State
State/State.cc
State/State.h

#Increment
Increment/Increment.cc
Increment/Increment.h

#Model
Model/Model.cc
Model/Model.h
Model/ModelFortran.h
Model/soca_model_mod.F90
Model/soca_model_interface_mod.F90
Model/soca_mom6.F90

#Linear Model
LinearModel/TlmId.cc
LinearModel/TlmId.h

#Covariance
Covariance/ErrorCovariance.cc
Covariance/ErrorCovariance.h
Covariance/soca_covariance_mod.F90
Covariance/soca_covariance_interface_mod.F90

#LocalizationMatrix

#Run
Run/Run.cc
Run/Run.h
Run/soca_setup.f90

#Utils
Utils/soca_utils.F90
Utils/random_vectors_mod.F90
Utils/random_f.cc
Utils/random_f.h
Utils/soca_bumpinterp2d_mod.F90
Utils/soca_omb_stats_mod.F90

#GetValuesTraj
GetValuesTraj/GetValuesTraj.cc
GetValuesTraj/GetValuesTraj.h
GetValuesTraj/soca_getvaltraj_mod.F90

#Variable transforms
Transforms/instantiateBalanceOpFactory.h
Transforms/util/soca_bkgerrutil_mod.F90

Transforms/Balance/Balance.cc
Transforms/Balance/Balance.h
Transforms/Balance/soca_balance_mod.F90
Transforms/Balance/soca_balance_interface_mod.F90
Transforms/Balance/soca_kst_mod.F90
Transforms/Balance/soca_ksshts_mod.F90

Transforms/BkgErr/BkgErr.cc
Transforms/BkgErr/BkgErr.h
Transforms/BkgErr/BkgErrFortran.h
Transforms/BkgErr/soca_bkgerr_mod.F90
Transforms/BkgErr/soca_bkgerr_interface_mod.F90

Transforms/BkgErrGodas/BkgErrGodas.cc
Transforms/BkgErrGodas/BkgErrGodas.h
Transforms/BkgErrGodas/BkgErrGodasFortran.h
Transforms/BkgErrGodas/soca_bkgerrgodas_mod.F90
Transforms/BkgErrGodas/soca_bkgerrgodas_interface_mod.F90

Transforms/BkgErrFilt/BkgErrFilt.cc
Transforms/BkgErrFilt/BkgErrFilt.h
Transforms/BkgErrFilt/BkgErrFiltFortran.h
Transforms/BkgErrFilt/soca_bkgerrfilt_mod.F90
Transforms/BkgErrFilt/soca_bkgerrfilt_interface_mod.F90

Transforms/VertConv/VertConv.cc
Transforms/VertConv/VertConv.h
Transforms/VertConv/soca_vertconv_mod.F90
Transforms/VertConv/soca_vertconv_interface_mod.F90

instantiateObsFilterFactory.h
Traits.h
Fortran.h

ModelBias.h
ModelBiasIncrement.h
ModelBiasCovariance.h
)

ecbuild_add_library( TARGET   soca
                     SOURCES  ${soca_src_files}
                     LIBS     ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES}
                              eckit eckit_mpi fckit oops ioda ufo fms mom6
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${SOCA_LINKER_LANGUAGE}
                    )
