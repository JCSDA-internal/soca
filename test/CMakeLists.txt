# This file sets up the ctests for SOCA
#
# Each ctest will be run in its own working directory (to keep tests from
# breaking when run in parallel), this working directory is available at
# "./soca/test/test_workdir/<test_name>". Each test working directy has a link
# to the following directories:
#
# - data_static - input data that does NOT come from any other ctest
# - data_generated - contains the output from each ctest, and can be used as
#                    input for another ctest
# - data_output - a convenience link to "./data_generated/<test_name>", this is
#                 created to easily differentiate between input and output data
#                 for a given ctest yaml


set( soca_test_input
  testinput/3dhyb_diffusion.yml
  testinput/3dhyb.yml
  testinput/3dhybfgat_pseudo.yml
  testinput/3dvar_moreice.yml
  testinput/3dvar_lowres.yml
  testinput/3dvar_nicas.yml
  testinput/3dvar_single_ob.yml
  testinput/3dvar_zero_ob.yml
  testinput/3dvar.yml
  testinput/3dvarfgat_pseudo.yml
  testinput/4denvar.yml
  testinput/4dhybenvar.yml
  testinput/4dvar_identity.yml
  testinput/4dvar_htlm.yml
  testinput/addincrement.yml
  testinput/convertincrement.yml
  testinput/convertstate_changevar.yml
  testinput/convertstate_soca2cice.yml
  testinput/convertstate.yml
  testinput/diffstates.yml
  testinput/dirac_diffusion.yml
  testinput/dirac_diffusion_lowres.yml
  testinput/dirac_mlbalance.yml
  testinput/dirac_soca_cor_nicas_scales.yml
  testinput/dirac_soca_cov.yml
  testinput/dirac_socahyb_cov.yml
  testinput/enshofx_1.yml
  testinput/enshofx_2.yml
  testinput/enshofx_3.yml
  testinput/enshofx.yml
  testinput/ensmeanandvariance.yml
  testinput/enspert.yml
  testinput/ensrecenter.yml
  testinput/forecast_identity.yml
  testinput/forecast_mom6_bgc.yml
  testinput/forecast_mom6_ens1.yml
  testinput/forecast_mom6_ens2.yml
  testinput/forecast_mom6_ens3.yml
  testinput/forecast_mom6.yml
  testinput/forecast_pseudo.yml
  testinput/gen_hybrid_linear_model_coeffs.yml
  testinput/geometry_iterator_2d.yml
  testinput/geometry_iterator_3d.yml
  testinput/geometry.yml
  testinput/getvalues.yml
  testinput/gridgen.yml
  testinput/hofx_3d.yml
  testinput/hofx_4d_pseudo.yml
  testinput/hofx_oasim_3d.yml
  testinput/hybridgain.yml
  testinput/increment.yml
  testinput/letkf_split_observer.yml
  testinput/letkf_split_solver.yml
  testinput/letkf.yml
  testinput/linearization_error.yml
  testinput/makeobs.yml
  testinput/model.yml
  testinput/obslocalization_vertical.yml
  testinput/obslocalization.yml
  testinput/obslocalizations.yml
  testinput/obsop_name_map.yml
  testinput/oceanSmoother.yml
  testinput/parameters_bump_cor_nicas_scales.yml
  testinput/parameters_bump_cor_nicas.yml
  testinput/parameters_bump_cov.yml
  testinput/parameters_bump_loc.yml
  testinput/parameters_diffusion.yml
  testinput/parameters_diffusion_lowres.yml
  testinput/parametric_stddev.yml
  testinput/setcorscales.yml
  testinput/sqrtvertloc.yml
  testinput/state.yml
  testinput/varchange_ana2model.yml
  testinput/varchange_balance_TSSSH.yml
  testinput/varchange_balance.yml
  testinput/varchange_bkgerrgodas.yml
  testinput/varchange_bkgerrsoca_stddev.yml
  testinput/varchange_bkgerrsoca.yml
)

set( soca_test_ref
  testref/3dhyb_diffusion.test
  testref/3dhyb.test
  testref/3dhybfgat_pseudo.test
  testref/3dvar_diffmlb.test
  testref/3dvar_lowres.test
  testref/3dvar_nicas.test
  testref/3dvar_single_ob.test
  testref/3dvar_zero_ob.test
  testref/3dvar.test
  testref/3dvar_moreice.test
  testref/3dvarfgat_pseudo.test
  testref/4denvar.test
  testref/4dhybenvar.test
  testref/4dvar_identity.test
  testref/4dvar_htlm.test
  testref/addincrement.test
  testref/convertincrement.test
  testref/convertstate_changevar.test
  testref/convertstate_soca2cice.test
  testref/convertstate.test
  testref/diffstates.test
  testref/dirac_diffusion.test
  testref/dirac_diffusion_lowres.test
  testref/dirac_soca_cor_nicas_scales.test
  testref/dirac_soca_cov.test
  testref/dirac_socahyb_cov.test
  testref/enshofx.test
  testref/ensmeanandvariance.test
  testref/enspert.test
  testref/ensrecenter.test
  testref/forecast_identity.test
  testref/forecast_pseudo.test
  testref/gen_hybrid_linear_model_coeffs.test
  testref/gridgen.test
  testref/hofx_3d.test
  testref/hofx_4d_pseudo.test
  testref/hofx_oasim_3d.test
  testref/hybridgain.test
  testref/letkf_split_observer.test
  testref/letkf_split_solver.test
  testref/letkf.test
  testref/linearization_error.test
  testref/makeobs.test
  testref/oceanSmoother.test
  testref/parameters_bump_cor_nicas_scales.test
  testref/parameters_bump_cor_nicas.test
  testref/parameters_bump_cov.test
  testref/parameters_bump_loc.test
  testref/parametric_stddev.test
  testref/setcorscales.test
  testref/sqrtvertloc.test
)

# Create directory for test input/output/reference and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testref)
foreach(FILENAME ${soca_test_input} ${soca_test_ref})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# create static data directory, and symlink all files.
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/Data
          ${CMAKE_CURRENT_BINARY_DIR}/data_static )
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data_generated)


# The files that need to be "installed" for outside repos, such as soca-science
#-------------------------------------------------------------------------------
set( soca_install_restarts
  Data/72x35x25/restarts/cice.res.nc
  Data/72x35x25/restarts/iced.2018-04-15-00000.nc
  Data/72x35x25/restarts/MOM.res.nc
)
install(FILES ${soca_install_restarts}
        DESTINATION ${INSTALL_DATA_DIR}/testdata/72x35x25/restarts/ )

set( soca_install_model_input
  Data/72x35x25/INPUT/forcing_daily.nc
  Data/72x35x25/INPUT/forcing_monthly.nc
  Data/72x35x25/INPUT/hycom1_25.nc
  Data/72x35x25/INPUT/layer_coord25.nc
  Data/72x35x25/INPUT/ocean_hgrid.nc
  Data/72x35x25/INPUT/ocean_topog.nc
)
install(FILES ${soca_install_model_input}
        DESTINATION ${INSTALL_DATA_DIR}/testdata/72x35x25/INPUT/ )

set( soca_install_data
  Data/rossrad.nc
  Data/godas_sst_bgerr.nc )
install(FILES ${soca_install_data}
        DESTINATION ${INSTALL_DATA_DIR}/testdata/ )

set( soca_install_mlbalance
  Data/mlbalance/antarctic.pt
  Data/mlbalance/normalization.antarctic.pt
  Data/mlbalance/arctic.pt
  Data/mlbalance/normalization.arctic.pt
)
install(FILES ${soca_install_mlbalance}
        DESTINATION ${INSTALL_DATA_DIR}/testdata/mlbalance/ )

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
# There are three types of tests:
#  1) soca executables (EXE must be given)
#  2) interface tests  (SRC must be given)
#  3) a test script    (SCRIPT must be given)
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
  set(singlevals NAME EXE SRC SCRIPT CFG MPI)
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

  # setup working directory
  # TODO, most of these are not needed for most tests, now that the MOM6
  # model has been mostly removed
  set( WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/test_workdir/${ARG_NAME})
  file( MAKE_DIRECTORY ${WORKDIR} )
  foreach(DIRNAME data_generated data_static testinput testoutput testref)
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_BINARY_DIR}/${DIRNAME}
      ${WORKDIR}/${DIRNAME} )
  endforeach()
  file(MAKE_DIRECTORY ${WORKDIR}/RESTART)
  foreach(FILENAME diag_table field_table INPUT)
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${WORKDIR}/data_static/workdir/${FILENAME}
      ${WORKDIR}/${FILENAME} )
  endforeach()

  # setup output directory
  file( MAKE_DIRECTORY ${WORKDIR}/data_generated/${ARG_NAME})
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${WORKDIR}/data_generated/${ARG_NAME}
    ${WORKDIR}/data_output
  )

  # Are we building a unit test / or running a soca executable?
  if ( ARG_SRC )
    # building a unit test
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      SOURCES executables/${ARG_SRC}
                      ARGS    ${CONFIG_FILE}
                      LIBS    soca
                      MPI     ${MPI}
                      ENVIRONMENT ${TRAPFPE_ENV}
                      WORKING_DIRECTORY ${WORKDIR}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
    set( EXE ./test_soca_${ARG_NAME})

  elseif ( ARG_EXE)
    # running a soca executable.
    set( EXE ${CMAKE_BINARY_DIR}/bin/${ARG_EXE})
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      MPI ${MPI}
                      COMMAND ${EXE}
                      WORKING_DIRECTORY ${WORKDIR}
                      ARGS ${CONFIG_FILE}
                      ENVIRONMENT
                              ${TRAPFPE_ENV}
                      DEPENDS ${ARG_EXE}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
  elseif ( ARG_SCRIPT )
    # running a test specific script, right now this is only for the mom6solo forecast
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${ARG_SCRIPT}
                      WORKING_DIRECTORY ${WORKDIR}
                      ARGS ${CONFIG_FILE}
                      DEPENDS mom6solo
                      ENVIRONMENT "BIN_DIR=${CMAKE_BINARY_DIR}/bin"
		                  "MPIEXE=${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPI}"
                      TEST_DEPENDS ${ARG_TEST_DEPENDS}
                      )
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



#---------------------------------------------------------------------------------------------------
# Gridgen, which pretty much every other test depends on. It's special.
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME gridgen
               NOTRAPFPE
               EXE  soca_gridgen.x )



#---------------------------------------------------------------------------------------------------
# Low level class interface tests
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME geometry
               SRC  TestGeometry.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME geometry_iterator_2d
               SRC  TestGeometryIterator.cc
               TEST_DEPENDS test_soca_gridgen )

# 3D geometry iterator is currently not working (or needed?)
# soca_add_test( NAME geometry_iterator_3d
#                SRC  TestGeometryIterator.cc
#                TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME state
               SRC  TestState.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME increment
               SRC  TestIncrement.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME modelaux
               SRC  TestModelAuxControl.cc
               CFG  model.yml
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME getvalues
               SRC  TestGetValues.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME oceanSmoother
               SRC  TestOceanSmoother.cc
               TEST_DEPENDS test_soca_gridgen )

# Variable changes
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME varchange_ana2model
               SRC  TestNLVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_balance
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_balance_TSSSH
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrsoca
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrsoca_stddev
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME varchange_bkgerrgodas
               SRC  TestVariableChange.cc
               TEST_DEPENDS test_soca_gridgen )

# Observation localization methods
#---------------------------------------------------------------------------------------------------
soca_add_test( NAME obslocalization
               SRC  TestObsLocalization.cc
               TEST_DEPENDS test_soca_gridgen )

## vertical localization is not working right now. h was removed from the
## geometry (it is part of the state), so we need to rethink how the
## geometryiterator works in the vertical
# soca_add_test( NAME obslocalization_vertical
#                SRC  TestObsLocalization.cc )

soca_add_test( NAME obslocalizations
               SRC  TestObsLocalizations.cc )



#---------------------------------------------------------------------------------------------------
# The following tests use the main SOCA applications
#  TEST output stream is compared against existing testoutput/ files
#---------------------------------------------------------------------------------------------------

# model forecasts
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME forecast_identity
               EXE  soca_forecast.x
               TEST_DEPENDS test_soca_gridgen)

# Test mom6solo forecast
soca_add_test( NAME forecast_mom6
               SCRIPT mom6solo.py )

# # Test mom6_bgc forecast (optional)
# if ( ENABLE_OCEAN_BGC )
# soca_add_test( NAME forecast_mom6_bgc
#                EXE  soca_forecast.x
#                NOTRAPFPE
#                TEST_DEPENDS test_soca_gridgen )
# endif()

# Test PseudoModel forecast
soca_add_test( NAME forecast_pseudo
               EXE  soca_forecast.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_forecast_mom6 )

# create ensemble forecast from mom6solo
soca_add_test( NAME forecast_mom6_ens1
               SCRIPT mom6solo.py )

soca_add_test( NAME forecast_mom6_ens2
               SCRIPT mom6solo.py )

soca_add_test( NAME forecast_mom6_ens3
               SCRIPT mom6solo.py )


# background error
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME setcorscales
               EXE  soca_setcorscales.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen)

soca_add_test( NAME parameters_bump_cor_nicas
               NOTRAPFPE
               EXE  soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen)

soca_add_test( NAME parameters_bump_cor_nicas_scales
               NOTRAPFPE
               EXE  soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_setcorscales)

soca_add_test( NAME parameters_bump_loc
               EXE  soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME parameters_bump_cov
               EXE  soca_error_covariance_toolbox.x
               NOTRAPFPE )

soca_add_test( NAME parameters_diffusion
               EXE  soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_setcorscales)

soca_add_test( NAME parameters_diffusion_lowres
               EXE  soca_error_covariance_toolbox.x)


# Misc applications
#---------------------------------------------------------------------------------------------------

# ensemble generation
soca_add_test( NAME enspert
               EXE  soca_enspert.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion )

# Remapping MOM6 (horiz+vertical intterpolation)
soca_add_test( NAME convertstate
               EXE  soca_convertstate.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Apply a nonlinear change of variable to an ensemble of states
soca_add_test( NAME convertstate_changevar
               EXE  soca_convertstate.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME ensmeanandvariance
               EXE  soca_ensmeanandvariance.x
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME parametric_stddev
               EXE  soca_convertincrement.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME ensrecenter
               EXE  soca_ensrecenter.x
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME hybridgain
               EXE  soca_hybridgain.x
               TEST_DEPENDS test_soca_gridgen )

# Generate sqrt of vertical localization
soca_add_test( NAME sqrtvertloc
               EXE  soca_sqrtvertloc.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion )

soca_add_test( NAME diffstates
               EXE  soca_diffstates.x
               TEST_DEPENDS test_soca_gridgen )

# Compute HTLM coefficients
soca_add_test( NAME gen_hybrid_linear_model_coeffs
               EXE  soca_gen_hybrid_linear_model_coeffs.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_forecast_mom6
                            test_soca_forecast_mom6_ens1
                            test_soca_forecast_mom6_ens2
                            test_soca_forecast_mom6_ens3 )

soca_add_test( NAME linearization_error
               EXE  soca_linearization_error.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_forecast_mom6
                            test_soca_forecast_mom6_ens1
                            test_soca_forecast_mom6_ens2
                            test_soca_forecast_mom6_ens3 )

# Dirac tests
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME dirac_soca_cor_nicas_scales
               EXE  soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_bump_cor_nicas_scales)

soca_add_test( NAME dirac_soca_cov
               EXE  soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion  )

soca_add_test( NAME dirac_socahyb_cov
               EXE  soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_bump_loc
                            test_soca_parameters_diffusion )

# NOTE the following tests throws a FPE only on clang
soca_add_test( NAME dirac_diffusion
               EXE soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion )

soca_add_test( NAME dirac_diffusion_lowres
               EXE soca_error_covariance_toolbox.x
               TEST_DEPENDS test_soca_parameters_diffusion_lowres )


if( ${daml_FOUND} )
  soca_add_test( NAME dirac_mlbalance
                 EXE soca_error_covariance_toolbox.x)
endif()

# h(x) executables
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME makeobs
               EXE  soca_hofx.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME hofx_3d
               EXE  soca_hofx3d.x
               TEST_DEPENDS test_soca_gridgen )

if( ${oasim_FOUND} )
  soca_add_test( NAME hofx_oasim_3d
                 EXE  soca_hofx3d.x
                 TEST_DEPENDS test_soca_gridgen )
endif( ${oasim_FOUND} )

soca_add_test( NAME hofx_4d_pseudo
               EXE  soca_hofx.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_forecast_mom6 )

# NOTE, enshofx is disabled because MOM6 always writes out some diagnostic
# netcdf files (which can't be disabled). This is normally not a problem, except
# with enshofx, because multiple instances of MOM6 are started under separate
# MPI_COMMs, causing problems when they all try to write to the same file at the
# same time. We'll need to take a closer look at MOM6 output when we start
# dealing with in-core MOM6

# soca_add_test( NAME enshofx
#                EXE  soca_enshofx.x
#                MPI 3
#                NOTRAPFPE
#                TEST_DEPENDS test_soca_gridgen )


# variational methods
# (unless specifically specified, it is assumed that correlation uses
#  EXPLICIT_DIFFUSION, and localization uses BUMP_NICAS )
#---------------------------------------------------------------------------------------------------

# # Following test does not work yet. IODA crashes
# soca_add_test( NAME 3dvar_zero_ob
#                EXE  soca_var.x
#                NOTRAPFPE
#                TEST_DEPENDS test_soca_static_socaerror_init )

# # Following test is breaking on intel 2022, beyond 1 iteration solver produces NaNs
# soca_add_test( NAME 3dvar_single_ob
#                EXE  soca_var.x
#                NOTRAPFPE
#                TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME 3dvar
               EXE  soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion )

soca_add_test( NAME 3dvar_moreice
               EXE  soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion )

soca_add_test( NAME 3dvar_nicas
               EXE  soca_var.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_parameters_bump_cor_nicas )


if( ${daml_FOUND} )
  soca_add_test( NAME 3dvar_diffmlb
                 EXE  soca_var.x
                 TEST_DEPENDS test_soca_gridgen
                              test_soca_parameters_diffusion )
endif()

soca_add_test( NAME 3dvar_lowres
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion_lowres )

soca_add_test( NAME 3dvarfgat_pseudo
               EXE  soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_forecast_mom6
                            test_socs_parameters_diffusion )

soca_add_test( NAME 3dhyb
               EXE  soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_bump_loc
                            test_socs_parameters_diffusion )

soca_add_test( NAME 3dhyb_diffusion
               EXE  soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_diffusion)

soca_add_test( NAME 3dhybfgat_pseudo
               EXE  soca_var.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_bump_loc
                            test_socs_parameters_diffusion )

soca_add_test( NAME 4denvar
               EXE soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_forecast_mom6_ens1
                            test_soca_forecast_mom6_ens2
                            test_soca_forecast_mom6_ens3 )

soca_add_test( NAME 4dhybenvar
               EXE soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_socs_parameters_diffusion
                            test_soca_forecast_mom6
                            test_soca_forecast_mom6_ens1
                            test_soca_forecast_mom6_ens2
                            test_soca_forecast_mom6_ens3 )

soca_add_test( NAME 4dvar_identity
               EXE soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_parameters_bump_loc
                            test_socs_parameters_diffusion
                            test_soca_forecast_mom6 )

soca_add_test( NAME 4dvar_htlm
               EXE soca_var.x
               TEST_DEPENDS test_soca_gridgen
                            test_socs_parameters_diffusion
                            test_soca_forecast_mom6
                            test_soca_forecast_mom6_ens1
                            test_soca_forecast_mom6_ens2
                            test_soca_forecast_mom6_ens3
                            test_soca_gen_hybrid_linear_model_coeffs )
# LETKF
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME letkf
               EXE  soca_letkf.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME letkf_split_observer
               EXE  soca_letkf.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME letkf_split_solver
               EXE  soca_letkf.x
              #  NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_letkf_split_observer )


# other executables that depends on increments from 3dvar_diffusion run
#---------------------------------------------------------------------------------------------------

soca_add_test( NAME addincrement
               EXE  soca_addincrement.x
               TEST_DEPENDS test_soca_gridgen
                            test_soca_3dvar )

# Remapping MOM6 (horiz+vertical interpolation)
soca_add_test( NAME convertincrement
               EXE  soca_convertincrement.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen
                            test_soca_3dvar)


if( ${icepack_FOUND} )
  # Apply the Soca2Cice nonlinear change of variable
  # TODO (G): check test dependency
  file(COPY ${CMAKE_CURRENT_BINARY_DIR}/data_static/72x35x25/restarts/iced.2018-04-15-00000.nc
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data_generated/convertstate_soca2cice/ )
     soca_add_test( NAME convertstate_soca2cice
                    EXE  soca_convertstate.x
                    MPI 2
                    NOTRAPFPE
                    TEST_DEPENDS test_soca_gridgen
                                 test_soca_3dvar )
endif()
