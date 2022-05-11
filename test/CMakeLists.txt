set( soca_test_input
  testinput/3dhyb.yml
  testinput/3dhybfgat.yml
  testinput/3dvarbump.yml
  testinput/3dvar_godas.yml
  testinput/3dvar_single_ob.yml
  testinput/3dvar_soca.yml
  testinput/3dvar_zero_ob.yml
  testinput/3dvarlowres_soca.yml
  testinput/3dvarfgat.yml
  testinput/3dvarfgat_pseudo.yml
  testinput/addincrement.yml
  testinput/balance_mask.yml
  testinput/checkpointmodel.yml
  testinput/convertincrement.yml
  testinput/convertstate.yml
  testinput/convertstate_changevar.yml
  testinput/diffstates.yml
  testinput/dirac_horizfilt.yml
  testinput/dirac_soca_cov.yml
  testinput/dirac_soca_mask.yml
  testinput/dirac_soca_nomask.yml
  testinput/dirac_socahyb_cov.yml
  testinput/dirac_soca_cor_nicas_scales.yml
  testinput/enshofx_1.yml
  testinput/enshofx_2.yml
  testinput/enshofx_3.yml
  testinput/enshofx.yml
  testinput/enspert.yml
  testinput/ensrecenter.yml
  testinput/ensvariance.yml
  testinput/forecast_identity.yml
  testinput/forecast_mom6.yml
  testinput/forecast_mom6_bgc.yml
  testinput/forecast_pseudo.yml
  testinput/geometry.yml
  testinput/geometry_iterator_2d.yml
  testinput/geometry_iterator_3d.yml
  testinput/gridgen.yml
  testinput/gridgen_small.yml
  testinput/hofx_3d.yml
  testinput/hofx_4d.yml
  testinput/hofx_4d_pseudo.yml
  testinput/hybridgain.yml
  testinput/increment.yml
  testinput/letkf_split_observer.yml
  testinput/letkf_split_solver.yml
  testinput/letkf.yml
  testinput/linearmodel.yml
  testinput/makeobs.yml
  testinput/model.yml
  testinput/model_ufsm6c6.yml
  testinput/obslocalization.yml
  testinput/obslocalization_vertical.yml
  testinput/obslocalizations.yml
  testinput/parameters_bump_cor_nicas.yml
  testinput/parameters_bump_cor_nicas_scales.yml
  testinput/parameters_bump_cov_lct.yml
  testinput/parameters_bump_cov_nicas.yml
  testinput/parameters_bump_loc.yml
  testinput/sqrtvertloc.yml
  testinput/setcorscales.yml
  testinput/state.yml
  testinput/static_socaerror_init.yml
  testinput/static_socaerrorlowres_init.yml
  testinput/static_socaerror_test.yml
  testinput/varchange_ana2model.yml
  testinput/varchange_balance.yml
  testinput/varchange_balance_TSSSH.yml
  testinput/varchange_bkgerrfilt.yml
  testinput/varchange_bkgerrgodas.yml
  testinput/varchange_bkgerrsoca.yml
  testinput/varchange_horizfilt.yml
  testinput/varchange_vertconv.yml
)

set( soca_test_ref
  testref/3dhyb.test
  testref/3dhybfgat.test
  testref/3dvarbump.test
  testref/3dvar_godas.test
  testref/3dvar_single_ob.test
  testref/3dvar_soca.test
  testref/3dvar_zero_ob.test
  testref/3dvarlowres_soca.test
  testref/3dvarfgat.test
  testref/3dvarfgat_pseudo.test
  testref/addincrement.test
  testref/balance_mask.test
  testref/checkpointmodel.test
  testref/convertincrement.test
  testref/convertstate.test
  testref/convertstate_changevar.test
  testref/diffstates.test
  testref/dirac_horizfilt.test
  testref/dirac_soca_cov.test
  testref/dirac_soca_mask.test
  testref/dirac_soca_nomask.test
  testref/dirac_socahyb_cov.test
  testref/dirac_soca_cor_nicas_scales.test
  testref/enshofx.test
  testref/enspert.test
  testref/ensrecenter.test
  testref/ensvariance.test
  testref/forecast_identity.test
  testref/forecast_mom6.test
  testref/forecast_mom6_bgc.test
  testref/forecast_pseudo.test
  testref/gridgen.test
  testref/gridgen_small.test
  testref/hofx_3d.test
  testref/hofx_4d.test
  testref/hofx_4d_pseudo.test
  testref/hybridgain.test
  testref/letkf_split_observer.test
  testref/letkf_split_solver.test
  testref/letkf.test
  testref/makeobs.test
  testref/parameters_bump_cor_nicas.test
  testref/parameters_bump_cor_nicas_scales.test
  testref/parameters_bump_cov_lct.test
  testref/setcorscales.test
  testref/sqrtvertloc.test
  testref/static_socaerror_init.test
  testref/static_socaerrorlowres_init.test
)

set( soca_data_files
  Data/godas_sst_bgerr.nc
  Data/rossrad.dat
  Data/fields_metadata.yml
)

set( soca_model_param
  Data/72x35x25/diag_table
  Data/72x35x25/field_table
  Data/72x35x25/ice.bkgerror.nc
  Data/72x35x25/ocn.bkgerror.nc
  Data/72x35x25/wav.bkgerror.nc
  Data/72x35x25/MOM_input
  Data/72x35x25/MOM_input_small
  Data/72x35x25/MOM_override_bgc
)

set( soca_model_namelist
  Data/72x35x25/input.nml
  Data/72x35x25/input_bgc.nml
  Data/72x35x25/input_small.nml
)

set( soca_model_restarts
  Data/72x35x25/INPUT/forcing_daily.nc
  Data/72x35x25/INPUT/forcing_monthly.nc
  Data/72x35x25/INPUT/hycom1_25.nc
  Data/72x35x25/INPUT/cice.res.nc
  Data/72x35x25/INPUT/layer_coord25.nc
  Data/72x35x25/INPUT/MOM.res.nc
  Data/72x35x25/INPUT/sfc.res.nc
  Data/72x35x25/INPUT/wav.res.nc
  Data/72x35x25/INPUT/ocean_hgrid.nc
  Data/72x35x25/INPUT/ocean_topog.nc
  Data/72x35x25/INPUT/ocean_mosaic.nc
  Data/72x35x25/INPUT/grid_spec.nc
)

set( soca_obs
  Data/Obs/adt.nc
  Data/Obs/biomass_p.nc
  Data/Obs/chl.nc
  Data/Obs/gmi_gpm_geoval.nc
  Data/Obs/gmi_gpm_obs.nc
  Data/Obs/icec.nc
  Data/Obs/icefb.nc
  Data/Obs/prof.nc
  Data/Obs/single_ob.nc
  Data/Obs/sss.nc
  Data/Obs/sst.nc
  Data/Obs/swh.nc
  Data/Obs/uocn_surface.nc
  Data/Obs/vocn_surface.nc
  Data/Obs/zero_ob.nc
)

# Create Data directory for test input/output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
foreach(FILENAME ${soca_test_input} ${soca_test_ref})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for model input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/INPUT)
foreach(FILENAME ${soca_model_restarts})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/INPUT/${filename} )
endforeach(FILENAME)
install(FILES ${soca_model_restarts}
        DESTINATION ${INSTALL_DATA_DIR}/testdata/72x35x25/INPUT/ )

# MOM's resource files
foreach(FILENAME ${soca_model_param} ${soca_data_files})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)
install(FILES ${soca_data_files}
        DESTINATION ${INSTALL_DATA_DIR}/testdata )


# FMS input.nml that may be overwritten during the testing
foreach(FILENAME ${soca_model_namelist})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/inputnml/${filename} )
endforeach(FILENAME)

# MOM's output directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RESTART)

ecbuild_add_resources( TARGET   soca_test_scripts
                       SOURCES_PACK
                       ${soca_test_input}
                       )
# Link to SOCA obs
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
foreach(FILENAME ${soca_obs})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )
endforeach(FILENAME)

# Create directory for BUMP output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump_lowres)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump-kmask)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump_scales)

#####################################################################

# Configurable parameters for the subsequent tests

set( SOCA_TESTS_MAX_MPI 2 CACHE STRING
  "Maximum number of PEs to use for SOCA's MPI based tests.\
  NOTE: should leave this at 2 to ensure all tests pass.")

set( SOCA_TESTS_FORC_TRAPFPE OFF CACHE BOOL
  "If true, the specified per-test trapping on/off flag is ignored\
  and is forced to ON for all tests" )

set( SOCA_TESTS_VALGRIND OFF CACHE BOOL
  "If true, some tests are run under valgrind")


#-------------------------------------------------------------------------------
# The following is a wrapper to simplify the generation of tests.
# There are two types of tests:
#  1) soca executables (EXE must be given)
#  2) interface tests  (SRC must be given)
# Soca executables by default will also run a second step where the output test
# log stream results are compared with reference values. The interface tests
# do not do this.
#
# Arguments:
#  NAME      - the name of the test (after test_soca is prepended), as well
#              as the base name of the input confg file
#  EXE       - name of executable to be used.
#  SRC       - The name of the unit test source file to compile
#  CFG       - The name of the yaml file, if the default testinput/${NAME}.yaml
#              is not to be used
#  MPI       - The number of MPI PEs to use. If not specified, ${SOCA_TESTS_MAX_MPI}
#              will be used
# TEST_DEPENDS - list of tests this test depends on
#-------------------------------------------------------------------------------
function(soca_add_test)
  # parse the passed arguments
  set(prefix     ARG)
  set(novals     NOTRAPFPE NOVALGRIND)
  set(singlevals NAME EXE SRC CFG MPI)
  set(multivals  TEST_DEPENDS)
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

  # determine if floating point error trapping should be set
  if ( ARG_NOTRAPFPE AND NOT SOCA_TESTS_FORC_TRAPFPE)
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=0")
  else()
    set ( TRAPFPE_ENV "OOPS_TRAPFPE=1")
  endif()

  if ( NOT SOCA_TESTS_VALGRIND )
    set ( ARG_NOVALGRIND TRUE)
  endif()

  # determine the default config file name
  if ( ARG_CFG )
    set ( CONFIG_FILE testinput/${ARG_CFG} )
  else()
    set ( CONFIG_FILE testinput/${ARG_NAME}.yml )
  endif()

  # MPI PEs
  set( MPI ${SOCA_TESTS_MAX_MPI} )
  if ( ARG_MPI )
    set( MPI ${ARG_MPI})
  endif()

  # Are we building a unit test / or running a soca executable?
  if ( ARG_SRC )
    # building a unit test
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      SOURCES executables/${ARG_SRC}
                      ARGS    ${CONFIG_FILE}
                      LIBS    soca
                      MPI     ${MPI}
                      ENVIRONMENT ${TRAPFPE_ENV}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
    set( EXE ./test_soca_${ARG_NAME})

  else()
    # running a soca executable.
    set( EXE ${CMAKE_BINARY_DIR}/bin/${ARG_EXE})
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      MPI ${MPI}
                      COMMAND ${EXE}
                      ARGS ${CONFIG_FILE}
                      ENVIRONMENT
                              ${TRAPFPE_ENV}
                      DEPENDS ${ARG_EXE}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})


  endif()

  # optional valgrind tests
  if ( NOT ARG_NOVALGRIND )
    ecbuild_add_test( TARGET test_soca_${ARG_NAME}_valgrind
                      MPI ${MPI}
                      COMMAND valgrind
                      ARGS -q --error-exitcode=42
                           --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/valgrind.sup
                           --gen-suppressions=all
                           ${EXE} ${CONFIG_FILE}
                      TEST_DEPENDS test_soca_${ARG_NAME}
                    )
  endif()
endfunction()



#-------------------------------------------------------------------------------
# Tests that create data other tests will use
#-------------------------------------------------------------------------------

# Create subsampled netcdf files
add_test ( test_soca_subsample_netcdf
           bash
           ${CMAKE_CURRENT_SOURCE_DIR}/Data/subsample_mom6files.sh )

# Test of grid generation and create subsampled grid
soca_add_test( NAME gridgen_small
               EXE  soca_gridgen.x
               TEST_DEPENDS test_soca_subsample_netcdf )

# Test of grid generation and create regular grid
soca_add_test( NAME gridgen
               EXE  soca_gridgen.x )

# Remapping MOM6 (horiz+vertical intterpolation)
soca_add_test( NAME convertstate
               EXE  soca_convertstate.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Apply a nonlinear change of variable to an ensemble of states
soca_add_test( NAME convertstate_changevar
               EXE  soca_convertstate.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Test of read/write restart for MOM6&ICE, model advance = Id
soca_add_test( NAME forecast_identity
               EXE  soca_forecast.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen)

# Test mom6 forecast
soca_add_test( NAME forecast_mom6
               EXE  soca_forecast.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Test mom6_bgc forecast (optional)
if ( ENABLE_OCEAN_BGC )
soca_add_test( NAME forecast_mom6_bgc
               EXE  soca_forecast.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )
endif()

# Test PseudoModel forecast
soca_add_test( NAME forecast_pseudo
               EXE  soca_forecast.x
               TEST_DEPENDS test_soca_forecast_mom6 )

# Initialize static error covariance model defined in soca, used for 3DVAR
soca_add_test( NAME static_socaerror_init
               EXE  soca_staticbinit.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME static_socaerrorlowres_init
               EXE  soca_staticbinit.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Use the soca_dirac.x application to create a mask for the dynamic height Jacobians
soca_add_test( NAME balance_mask
               EXE  soca_dirac.x
               TEST_DEPENDS test_soca_gridgen )

# Create dc/dt file for ocean/ice balance and convert dirac file into a mask
add_test ( test_soca_create_kmask
           python3
           ${CMAKE_CURRENT_SOURCE_DIR}/Data/create_kmask.py
           TEST_DEPENDS test_static_socaerror_init
                        test_balance_mask )

# ensemble generation
soca_add_test( NAME enspert
               EXE  soca_enspert.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_static_socaerror_init
                            test_soca_balance_mask
                            test_soca_create_kmask)

# Generate sqrt of vertical localization
soca_add_test( NAME sqrtvertloc
               EXE  soca_sqrtvertloc.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_static_socaerror_init
                            test_soca_balance_mask
                            test_soca_create_kmask)

# TODO Test below can take a few hours to complete.
# Revive once we have tiered testing?
#soca_add_test( NAME covar_lengthscale
#               EXE  soca_staticbinit.x
#               TEST_DEPENDS test_soca_gridgen)

# Generate synthetic observations
soca_add_test( NAME makeobs
               EXE  soca_hofx.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Tests of class interfaces
#-------------------------------------------------------------------------------
soca_add_test( NAME geometry
               SRC  TestGeometry.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME geometry_iterator_2d
               SRC  TestGeometryIterator.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME geometry_iterator_3d
               SRC  TestGeometryIterator.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME state
               SRC  TestState.cc
               TEST_DEPENDS test_soca_gridgen )


# Test MOM6-solo
soca_add_test( NAME model
               SRC  TestModel.cc
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME modelaux
               SRC  TestModelAuxControl.cc
               CFG  model.yml
               TEST_DEPENDS test_soca_gridgen )

# Test the UFS-mom6-cice6
if (BUILD_UFSM6C6)
  soca_add_test( NAME model_ufsm6c6
                 SRC  TestModel.cc
                 NOTRAPFPE
                 TEST_DEPENDS test_soca_gridgen )

  soca_add_test( NAME modelaux_ufsm6c6
                 SRC  TestModelAuxControl.cc
                 CFG  model_ufsm6c6.yml
                 TEST_DEPENDS test_soca_gridgen )
endif ()

soca_add_test( NAME increment
               SRC  TestIncrement.cc
               TEST_DEPENDS test_soca_gridgen )


soca_add_test( NAME errorcovariance
               SRC  TestErrorCovariance.cc
               CFG  static_socaerror_test.yml
               TEST_DEPENDS test_soca_gridgen
                            test_soca_static_socaerror_init )

soca_add_test( NAME linearmodel
               SRC  TestLinearModel.cc
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_static_socaerror_init )

soca_add_test( NAME varchange_ana2model
               SRC  TestNLVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_balance
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen
                            test_soca_create_dcdt )

soca_add_test( NAME varchange_balance_TSSSH
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrfilt
               SRC TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_horizfilt
               SRC TestVariableChange.cc
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrsoca
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrgodas
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_vertconv
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME obslocalization
               SRC  TestObsLocalization.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME obslocalization_vertical
               SRC  TestObsLocalization.cc
               TEST_DEPENDS test_soca_gridgen_small )

soca_add_test( NAME obslocalizations
               SRC  TestObsLocalizations.cc
               TEST_DEPENDS test_soca_gridgen_small )


#-------------------------------------------------------------------------------
# Tests using the main SOCA applications
#  TEST output stream is compared against existing testoutput/ files
#-------------------------------------------------------------------------------
soca_add_test( NAME ensvariance
               EXE  soca_ensvariance.x
               TEST_DEPENDS test_soca_enspert )

soca_add_test( NAME parameters_bump_loc
               EXE  soca_error_covariance_training.x
               TEST_DEPENDS test_soca_enspert )

soca_add_test( NAME ensrecenter
               EXE  soca_ensrecenter.x
               TEST_DEPENDS test_soca_enspert )

soca_add_test( NAME hybridgain
               EXE  soca_hybridgain.x
               TEST_DEPENDS test_soca_enspert )

# TODO the next couple of tests are a little slow.
# Revive them once we have tiered testing?
# soca_add_test( NAME parameters_bump_cov_lct
#                EXE  soca_error_covariance_training.x
#                TEST_DEPENDS test_soca_enspert )
#
# soca_add_test( NAME parameters_bump_cov_nicas
#                EXE  soca_error_covariance_training.x
#                TEST_DEPENDS test_soca_parameters_bump_cov_lct )

# soca_add_test( NAME dirac_bump_cov
#                EXE  soca_dirac.x )

soca_add_test( NAME setcorscales
               EXE  soca_setcorscales.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen)

soca_add_test( NAME parameters_bump_cor_nicas
               NOTRAPFPE
               EXE  soca_error_covariance_training.x
               TEST_DEPENDS test_soca_gridgen_small
                            test_soca_convertstate)

soca_add_test( NAME parameters_bump_cor_nicas_scales
               NOTRAPFPE
               EXE  soca_error_covariance_training.x
               TEST_DEPENDS test_soca_setcorscales)


# Dirac tests

soca_add_test( NAME dirac_soca_cor_nicas_scales
               EXE  soca_dirac.x
               TEST_DEPENDS parameters_bump_cor_nicas_scales)

soca_add_test( NAME dirac_soca_cov
               EXE  soca_dirac.x
               TEST_DEPENDS test_soca_static_socaerror_init)

soca_add_test( NAME dirac_socahyb_cov
               EXE  soca_dirac.x
               TEST_DEPENDS test_soca_parameters_bump_loc)

soca_add_test( NAME dirac_horizfilt
               EXE  soca_dirac.x
               TEST_DEPENDS test_soca_static_socaerror_init)

soca_add_test( NAME dirac_soca_mask
               EXE  soca_dirac.x)

soca_add_test( NAME dirac_soca_nomask
               EXE  soca_dirac.x)

#ecbuild_add_test( TARGET test_soca_obserrorcov
#                  SOURCES executables/TestObsErrorCovariance.cc
#                  ARGS "testinput/obserrorcovariance_test.yml"
#                  LIBS    soca )
#

# h(x) executables
soca_add_test( NAME hofx_3d
               EXE  soca_hofx3d.x
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME hofx_4d
               EXE  soca_hofx.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME hofx_4d_pseudo
               EXE  soca_hofx.x
               TEST_DEPENDS test_soca_forecast_mom6 )

# Re-enable the comparison once that's fixed
soca_add_test( NAME enshofx
               EXE  soca_enshofx.x
               MPI 3
               NOTRAPFPE
               TEST_DEPENDS test_soca_enspert )


# variational methods
soca_add_test( NAME 3dvar_soca
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dvar_single_ob
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

# # Following test does not work yet. IODA crashes
# soca_add_test( NAME 3dvar_zero_ob
#                EXE  soca_var.x
#                TOL  4e-9 0
#                NOTRAPFPE
#                TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dvarbump
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_parameters_bump_cor_nicas )

soca_add_test( NAME 3dvar_godas
               EXE  soca_var.x
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME addincrement
               EXE  soca_addincrement.x
               TEST_DEPENDS 3dvar_godas )

# Remapping MOM6 (horiz+vertical intterpolation)
soca_add_test( NAME convertincrement
               EXE  soca_convertincrement.x
               TOL  2e-5 0
               NOTRAPFPE
               TEST_DEPENDS 3dvar_godas )

soca_add_test( NAME 3dvarlowres_soca
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerrorsmall_init )

soca_add_test( NAME diffstates
               EXE  soca_diffstates.x
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME 3dvarfgat
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dvarfgat_pseudo
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dhyb
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init
                            test_soca_parameters_bump_loc )

soca_add_test( NAME 3dhybfgat
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init
                            test_soca_parameters_bump_loc )

# soca_add_test( NAME 4dhyb
#                EXE  soca_var.x
#                TEST_DEPENDS test_soca_static_socaerror_init
#                             test_soca_parameters_bump_loc )

soca_add_test( NAME letkf
               EXE  soca_letkf.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_enspert )

soca_add_test( NAME letkf_split_observer
               EXE  soca_letkf.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_enspert )

ecbuild_add_test( TARGET letkf_split_observer_post
                  TYPE SCRIPT
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/letkf_obscat.sh
                  TEST_DEPENDS letkf_split_observer)

soca_add_test( NAME letkf_split_solver
               EXE  soca_letkf.x
               NOTRAPFPE
               TEST_DEPENDS letkf_split_observer_post )



# restart checkpointing
soca_add_test( NAME checkpointmodel
               EXE  soca_checkpoint_model.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_3dvar_godas)
