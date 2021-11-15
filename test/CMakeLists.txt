set( soca_test_input
  testinput/3dhyb.yml
  testinput/3dhybfgat.yml
  testinput/3dvarbump.yml
  testinput/3dvar_godas.yml
  testinput/3dvar_soca.yml
  testinput/3dvarlowres_soca.yml
  testinput/3dvarfgat.yml
  testinput/3dvarfgat_pseudo.yml
  testinput/addincrement.yml
  testinput/balance_mask.yml
  testinput/checkpointmodel.yml
  testinput/convertstate.yml
  testinput/convertstate_changevar.yml
  testinput/diffstates.yml
  testinput/dirac_horizfilt.yml
  testinput/dirac_soca_cov.yml
  testinput/dirac_socahyb_cov.yml
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
  testinput/geometry_iterator.yml
  testinput/geometryatm.yml
  testinput/getvalues.yml
  testinput/gridgen.yml
  testinput/gridgen_small.yml
  testinput/hofx_3d.yml
  testinput/hofx_3dcrtm.yml
  testinput/hofx_4d.yml
  testinput/hofx_4d_pseudo.yml
  testinput/hybridgain.yml
  testinput/increment.yml
  testinput/letkf_observer.yml
  testinput/letkf_solver.yml
  testinput/lineargetvalues.yml
  testinput/linearmodel.yml
  testinput/makeobs.yml
  testinput/model.yml
  testinput/obslocalization.yml
  testinput/parameters_bump_cor_nicas.yml
  testinput/parameters_bump_cov_lct.yml
  testinput/parameters_bump_cov_nicas.yml
  testinput/parameters_bump_loc.yml
  testinput/sqrtvertloc.yml
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
  testref/3dvar_soca.test
  testref/3dvarlowres_soca.test
  testref/3dvarfgat.test
  testref/3dvarfgat_pseudo.test
  testref/addincrement.test
  testref/balance_mask.test
  testref/checkpointmodel.test
  testref/convertstate.test
  testref/convertstate_changevar.test
  testref/diffstates.test
  testref/dirac_horizfilt.test
  testref/dirac_soca_cov.test
  testref/dirac_socahyb_cov.test
  testref/enshofx.test
  testref/enspert.test
  testref/ensrecenter.test
  testref/ensvariance.test
  testref/forecast_identity.test
  testref/forecast_mom6.test
  testref/forecast_mom6_bgc.test
  testref/forecast_pseudo.test
  testref/hofx_3d.test
  testref/hofx_3dcrtm.test
  testref/hofx_4d.test
  testref/hofx_4d_pseudo.test
  testref/hybridgain.test
  testref/letkf_observer.test
  testref/letkf_solver.test
  testref/makeobs.test
  testref/parameters_bump_cor_nicas.test
  testref/parameters_bump_cov_lct.test
  testref/sqrtvertloc.test
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
  Data/Obs/sss.nc
  Data/Obs/sst.nc
  Data/Obs/swh.nc
  Data/Obs/uocn_surface.nc
  Data/Obs/vocn_surface.nc
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


# Link CRTM coeff
if( ${crtm_FOUND} )

  # the coefficient files needed
  # TODO(someone) for now these are our own copy in git-lfs,
  #  because why download 627MB if I only need 7MB of files??
  #  Download these directly from aws instead when that is an option.
  list( APPEND crtm_coef
    Data/crtm/gmi_gpm.SpcCoeff.bin
    Data/crtm/gmi_gpm.TauCoeff.bin
    Data/crtm/CloudCoeff.bin
    Data/crtm/AerosolCoeff.bin
    Data/crtm/FASTEM6.MWwater.EmisCoeff.bin
  )

  # link the files
  foreach(FILENAME ${crtm_coef})
       get_filename_component(filename ${FILENAME} NAME )
       execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
             ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
             ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )
  endforeach(FILENAME)
endif()


# Create directory for BUMP output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump_lowres)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump-kmask)

# find compare.py provided by oops
find_file(COMPARE_BIN oops_compare.py HINTS ${oops_BINDIR})

#####################################################################

# Configurable parameters for the subsequent tests

set( SOCA_TESTS_MAX_MPI 2 CACHE STRING
  "Maximum number of PEs to use for SOCA's MPI based tests.\
  NOTE: should leave this at 2 to ensure all tests pass.")

set( SOCA_TESTS_DEFAULT_TOL "1.0e-12;0" CACHE STRING
  "Default test tolerances for compare.py based tests. (float_tol;int_tol) ")

set( SOCA_TESTS_FORC_DEFAULT_TOL OFF CACHE BOOL
  "If true, the specified test tolerances are ignored and\
  SOCA_TESTS_DEFAULT_TOL is used" )

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
#  NOCOMPARE - if flag is present, no compare test is performed
#  NAME      - the name of the test (after test_soca is prepended), as well
#              as the base name of the input confg file
#  EXE       - name of executable to be used.
#  SRC       - The name of the unit test source file to compile
#  CFG       - The name of the yaml file, if the default testinput/${NAME}.yaml
#              is not to be used
#  TOL       - The tolerances to use for the test "float_tolerance int_tolerance"
#              If none are given, ${SOCA_TESTS_DEFAULT_TOL} is used. Only used
#              if an EXE with no NOCOMPARE flag set.
#  MPI       - The number of MPI PEs to use. If not specified, ${SOCA_TESTS_MAX_MPI}
#              will be used
# TEST_DEPENDS - list of tests this test depends on
#-------------------------------------------------------------------------------
function(soca_add_test)
  # parse the passed arguments
  set(prefix     ARG)
  set(novals     NOCOMPARE NOTRAPFPE NOVALGRIND)
  set(singlevals NAME EXE SRC CFG MPI)
  set(multivals  TEST_DEPENDS TOL)
  cmake_parse_arguments(${prefix}
                        "${novals}" "${singlevals}" "${multivals}"
                        ${ARGN})

  # set the tolerances to defaults if none given, or if we want to override
  IF(NOT ARG_TOL OR SOCA_TESTS_FORC_DEFAULT_TOL)
    SET(ARG_TOL ${SOCA_TESTS_DEFAULT_TOL})
  ENDIF()

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
    # building a unit test, therfore also assume no compare step
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      SOURCES executables/${ARG_SRC}
                      ARGS    ${CONFIG_FILE}
                      LIBS    soca
                      MPI     ${MPI}
                      ENVIRONMENT ${TRAPFPE_ENV}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
    set( EXE ./test_soca_${ARG_NAME})

  else()
    LIST(GET ARG_TOL 0 TOL_F)
    LIST(GET ARG_TOL 1 TOL_I)

    # find the MPI command
    set(MPI_CMD "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI}")

    # running a soca executable.
    # This is run with the run wrapper, which will also optionally
    # run the comapre.py script afterwards
    ecbuild_add_test( TARGET  test_soca_${ARG_NAME}
                      TYPE    SCRIPT
                      COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/test_wrapper.sh"
                      ARGS    "${CMAKE_BINARY_DIR}/bin/${ARG_EXE}"
                              ${CONFIG_FILE}
                      ENVIRONMENT
                              ${TRAPFPE_ENV}
                              COMPARE_SCRIPT=${COMPARE_BIN}
                              COMPARE_TESTNAME=${ARG_NAME}
                              COMPARE_TOL_F=${TOL_F}
                              COMPARE_TOL_I=${TOL_I}
                              MPI_CMD=${MPI_CMD}
                              SKIP_COMPARE=${ARG_NOCOMPARE}
                      DEPENDS ${ARG_EXE}
                      TEST_DEPENDS ${ARG_TEST_DEPENDS})
    set( EXE ${CMAKE_BINARY_DIR}/bin/${ARG_EXE})

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
               NOCOMPARE
               TEST_DEPENDS test_soca_subsample_netcdf )

# Test of grid generation and create regular grid
soca_add_test( NAME gridgen
               EXE  soca_gridgen.x
               NOCOMPARE )

# Remapping MOM6 (horiz+vertical intterpolation)
soca_add_test( NAME convertstate
               EXE  soca_convertstate.x
               TOL  2e-5 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

# Apply a nonlinear change of variable to an ensemble of states
soca_add_test( NAME convertstate_changevar
               EXE  soca_convertstate.x
               TOL  2e-5 0
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
               NOCOMPARE
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME static_socaerrorlowres_init
               EXE  soca_staticbinit.x
               NOCOMPARE
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

soca_add_test( NAME geometry_iterator
               SRC  TestGeometryIterator.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME geometryatm
               SRC  TestGeometry.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME state
               SRC  TestState.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME getvalues
               SRC  TestGetValues.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME model
               SRC  TestModel.cc
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME modelaux
               SRC  TestModelAuxControl.cc
               CFG  model.yml
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME increment
               SRC  TestIncrement.cc
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME lineargetvalues
               SRC  TestLinearGetValues.cc
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
               TEST_DEPENDS test_soca_gridgen
               NOCOMPARE )


#-------------------------------------------------------------------------------
# Tests using the main SOCA applications
#  TEST output stream is compared against existing testoutput/ files
#-------------------------------------------------------------------------------
soca_add_test( NAME ensvariance
               EXE  soca_ensvariance.x
               TEST_DEPENDS test_soca_enspert )

soca_add_test( NAME parameters_bump_loc
               EXE  soca_parameters.x
               NOCOMPARE
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
#                EXE  soca_parameters.x
#                TEST_DEPENDS test_soca_enspert )
#
# soca_add_test( NAME parameters_bump_cov_nicas
#                EXE  soca_parameters.x
#                TEST_DEPENDS test_soca_parameters_bump_cov_lct )

# soca_add_test( NAME dirac_bump_cov
#                EXE  soca_dirac.x )

soca_add_test( NAME parameters_bump_cor_nicas
               NOTRAPFPE
               EXE  soca_parameters.x
               TEST_DEPENDS test_soca_gridgen_small
                            test_soca_convertstate)


# Dirac tests

soca_add_test( NAME dirac_soca_cov
               EXE  soca_dirac.x
               TEST_DEPENDS test_soca_static_socaerror_init)

soca_add_test( NAME dirac_socahyb_cov
               EXE  soca_dirac.x
               TEST_DEPENDS test_soca_parameters_bump_loc)

soca_add_test( NAME dirac_horizfilt
               EXE  soca_dirac.x
               TEST_DEPENDS test_soca_static_socaerror_init)

#ecbuild_add_test( TARGET test_soca_obserrorcov
#                  SOURCES executables/TestObsErrorCovariance.cc
#                  ARGS "testinput/obserrorcovariance_test.yml"
#                  LIBS    soca )
#

# h(x) executables
soca_add_test( NAME hofx_3d
               EXE  soca_hofx3d.x
               TEST_DEPENDS test_soca_gridgen )

# CRTM is suddenly throwing div by 0 errors, disable trapping
if( ${crtm_FOUND} )
 soca_add_test( NAME hofx_3dcrtm
                EXE  soca_hofx3d.x
                NOTRAPFPE
                TEST_DEPENDS test_soca_gridgen )
endif()

soca_add_test( NAME hofx_4d
               EXE  soca_hofx.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME hofx_4d_pseudo
               EXE  soca_hofx.x
               TEST_DEPENDS test_soca_forecast_mom6 )

# TODO enshofx currently breaks compare.py
# Re-enable the comparison once that's fixed
soca_add_test( NAME enshofx
               EXE  soca_enshofx.x
               MPI 3
               NOCOMPARE
               NOTRAPFPE
               TEST_DEPENDS test_soca_enspert )


# variational methods
soca_add_test( NAME 3dvar_soca
               EXE  soca_var.x
               TOL  4e-9 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dvarbump
               EXE  soca_var.x
               TOL  4e-9 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_parameters_bump_cor_nicas )

soca_add_test( NAME 3dvar_godas
               EXE  soca_var.x
               TOL  2e-6 0
               TEST_DEPENDS test_soca_static_socaerror_init )

# TODO: Remove NOCOMPARE when oops PR#585 is merged
soca_add_test( NAME addincrement
               EXE  soca_addincrement.x
               NOCOMPARE
               TEST_DEPENDS 3dvar_godas )

soca_add_test( NAME 3dvarlowres_soca
               EXE  soca_var.x
               TOL  0.1 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerrorsmall_init )

soca_add_test( NAME diffstates
               EXE  soca_diffstates.x
               TEST_DEPENDS test_soca_gridgen )

soca_add_test( NAME 3dvarfgat
               EXE  soca_var.x
               TOL  1e-6 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dvarfgat_pseudo
               EXE  soca_var.x
               TOL  1e-6 0
               NOTRAPFPE
               TEST_DEPENDS test_soca_static_socaerror_init )

soca_add_test( NAME 3dhyb
               EXE  soca_var.x
               TOL  2e-7 0
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

soca_add_test( NAME letkf_observer
               EXE  soca_letkf.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_enspert )

ecbuild_add_test( TARGET letkf_observer_post
                  TYPE SCRIPT
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/letkf_obscat.sh
                  TEST_DEPENDS letkf_observer)

soca_add_test( NAME letkf_solver
               EXE  soca_letkf.x
               NOTRAPFPE
               TEST_DEPENDS letkf_observer_post )

# restart checkpointing
soca_add_test( NAME checkpointmodel
               EXE  soca_checkpoint_model.x
               NOTRAPFPE
               TEST_DEPENDS test_soca_3dvar_godas)
