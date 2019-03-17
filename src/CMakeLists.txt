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

#Increment
./Increment/Increment.cc
./Increment/Increment.h

#Model
./Model/Model.cc
./Model/Model.h
./Model/soca_model.interface.F90
./Model/soca_model_mod.F90
./Model/soca_mom6.F90

#Linear Model
./LinearModel/TlmId.cc
./LinearModel/TlmId.h

#Covariance
./Covariance/ErrorCovariance.cc
./Covariance/ErrorCovariance.h
./Covariance/soca_covariance.interface.F90
./Covariance/soca_covariance_mod.F90

#LocalizationMatrix

#Run
./Run/Run.cc
./Run/Run.h
./Run/soca_setup.f90

#Utils
./Utils/soca_utils.F90
./Utils/random_vectors_mod.F90
./Utils/random_f.cc
./Utils/random_f.h
./Utils/soca_bumpinterp2d_mod.F90
./Utils/soca_omb_stats_mod.F90

#GetValuesTraj
./GetValuesTraj/GetValuesTraj.cc
./GetValuesTraj/GetValuesTraj.h
./GetValuesTraj/soca_getvaltraj_mod.F90

#Variable transforms
./Transforms/instantiateBalanceOpFactory.h
./Transforms/util/soca_bkgerrutil_mod.F90


./Transforms/Balance/Balance.cc
./Transforms/Balance/Balance.h
./Transforms/Balance/soca_balance.interface.F90
./Transforms/Balance/soca_balance_mod.F90
./Transforms/Balance/soca_kst_mod.F90
./Transforms/Balance/soca_ksshts_mod.F90

./Transforms/BkgErr/BkgErr.cc
./Transforms/BkgErr/BkgErr.h
./Transforms/BkgErrGodas/BkgErrFortran.h
./Transforms/BkgErr/soca_bkgerr.interface.F90
./Transforms/BkgErr/soca_bkgerr_mod.F90

./Transforms/BkgErrGodas/BkgErrGodas.cc
./Transforms/BkgErrGodas/BkgErrGodas.h
./Transforms/BkgErrGodas/BkgErrGodasFortran.h
./Transforms/BkgErrGodas/soca_bkgerrgodas.interface.F90
./Transforms/BkgErrGodas/soca_bkgerrgodas_mod.F90

./Transforms/BkgErrFilt/BkgErrFilt.cc
./Transforms/BkgErrFilt/BkgErrFilt.h
./Transforms/BkgErrFilt/BkgErrFiltFortran.h
./Transforms/BkgErrFilt/soca_bkgerrfilt.interface.F90
./Transforms/BkgErrFilt/soca_bkgerrfilt_mod.F90

./Transforms/VertConv/VertConv.cc
./Transforms/VertConv/VertConv.h
./Transforms/VertConv/soca_vertconv.interface.F90
./Transforms/VertConv/soca_vertconv_mod.F90

instantiateObsFilterFactory.h
Traits.h
Fortran.h

ModelBias.h
ModelBiasIncrement.h
ModelBiasCovariance.h
)

ecbuild_add_library( TARGET   soca
                     SOURCES  ${oops_soca_src_files}
                     LIBS     oops ioda ufo fms mom6 ${LAPACK_LIBRARIES} ${NETCDF_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

# add soca coding norms test
ecbuild_add_test( TARGET soca_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )
