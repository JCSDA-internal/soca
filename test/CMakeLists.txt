list( APPEND soca_test_input
  testinput/truth.yml
  testinput/forecast.yml 
  testinput/hofx_makeobs3d.yml
  testinput/state_test.yml
  testinput/model_test.yml
  testinput/linearmodel_test.yml
  testinput/geometry_test.yml
  testinput/increment_test.yml
  testinput/errorcovariance_test.yml
  testinput/obserrorcovariance_test.yml
  testinput/bump_static_b.yml  
  testinput/3dvar_test.yml
  testinput/3dvarfgat_test.yml
  testinput/localization_test.yml
  testinput/nicas_test.yml
  testinput/4denvar_nicas_test.yml
  #testinput/variablechange.yml
  testinput/ksshts.yml
  testinput/kst.yml
  testinput/ktc.yml
  testinput/bkgerr.yml
  testinput/vertconv.yml  
  #testinput/parametersbump.cov.yml
  testinput/static_scale_binit.yml
  testinput/static_binit.yml  
  compare.sh
  run.sh
)

list( APPEND soca_testoutput
  #testoutput/dirac.nicas.test
  testoutput/forecast.test
  testoutput/truth.test
  testoutput/hofx.test  
  testoutput/makeobs3d.test
  testoutput/3dvar.test  
)

list( APPEND soca_model_param
  Data/360x210x63/input.nml
  Data/360x210x63/MOM_input
  Data/360x210x63/MOM_override
  Data/360x210x63/Vertical_coordinate.nc
  Data/rossrad.dat
)

list( APPEND soca_model_restarts
  Data/360x210x63/INPUT/MOM.res.nc
  Data/360x210x63/INPUT/ice_model.res.nc
  Data/360x210x63/INPUT/ocean_hgrid.nc
  Data/360x210x63/INPUT/topog.nc
)

# Create Data directory for test input and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${soca_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create Data directory for reference output and symlink all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
foreach(FILENAME ${soca_testoutput})
    message(${FILENAME})
    message(${CMAKE_CURRENT_SOURCE_DIR})    
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

#Link to MOM's resource files
foreach(FILENAME ${soca_model_param})
     get_filename_component(filename ${FILENAME} NAME )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

ecbuild_add_resources( TARGET   soca_test_scripts
                       SOURCES_PACK
                       ${soca_test_input}
                     )
# Link to IODA Marine obs
file(GLOB ioda_obs "${ioda_SOURCE_DIR}/test/testinput/marine/*.nc*")
#list( APPEND soca_obs
#  ../../ufo/test/testinput/marineobs/Jason-3-2018-04-15.nc
#  ../../ufo/test/testinput/marineobs/seaice_obs-2018-04-15.nc
#  ../../ufo/test/testinput/marineobs/sea-ice-con-l2.nc
#  ../../ufo/test/testinput/marineobs/sst_obs-2018-04-15.nc
#  ../../ufo/test/testinput/marineobs/t0n156e_dy.nc
#  ../../ufo/test/testinput/marineobs/cryosat2-2018-04-15.nc
#)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)
foreach(FILENAME ${ioda_obs})
     get_filename_component(filename ${FILENAME} NAME )
     #configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
     #               ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} COPYONLY)
     #execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
     #      ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}    
     #      ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )
     execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/Data/${filename} )	   
endforeach(FILENAME)


#####################################################################

# Tests that create data other tests might use (should check what is really needed...)

ecbuild_add_test( TARGET test_soca_truth
                  TYPE SCRIPT
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_forecast.x testinput/truth.yml"
                  DEPENDS soca_forecast.x )

ecbuild_add_test( TARGET test_soca_forecast
                  TYPE SCRIPT
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_forecast.x testinput/forecast.yml"
                  DEPENDS soca_forecast.x )

ecbuild_add_test( TARGET test_soca_makeobs3d
                  TYPE SCRIPT
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_makeobs.x testinput/hofx_makeobs3d.yml"
                  DEPENDS soca_makeobs.x )

ecbuild_add_test( TARGET test_soca_covar_lengthscale
                  TYPE SCRIPT
		  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x testinput/static_scale_binit.yml"
                  DEPENDS soca_staticbinit.x )

ecbuild_add_test( TARGET test_soca_binit
                  TYPE SCRIPT
		  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_staticbinit.x testinput/static_binit.yml"
                  DEPENDS soca_staticbinit.x )


# Tests of class interfaces

ecbuild_add_test( TARGET  test_soca_geometry
                  BOOST
		  MPI     1
                  SOURCES executables/TestGeometry.cc
                  ARGS    "testinput/geometry_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_linearmodel
                  BOOST
		  MPI     1		  
                  SOURCES executables/TestLinearModel.cc
                  ARGS    "testinput/linearmodel_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_state
                  BOOST
                  SOURCES executables/TestState.cc
                  ARGS    "testinput/state_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_modelaux
                  BOOST
                  SOURCES executables/TestModelAuxControl.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_model
                  BOOST
		  MPI     1
                  SOURCES executables/TestModel.cc
                  ARGS    "testinput/model_test.yml"
                  LIBS    soca )
		  
ecbuild_add_test( TARGET  test_soca_increment
                  BOOST
		  MPI     1 		  
                  SOURCES executables/TestIncrement.cc
                  ARGS    "testinput/increment_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_errorcovariance
                  BOOST
                  SOURCES executables/TestErrorCovariance.cc
                  ARGS    "testinput/errorcovariance_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET  test_soca_localization
                  BOOST
                  SOURCES executables/TestLocalization.cc
                  ARGS    "testinput/localization_test.yml"
                  MPI     1
                  LIBS    soca )
		  
ecbuild_add_test( TARGET test_soca_obserrorcov
                  BOOST
                  SOURCES executables/TestObsErrorCovariance.cc
                  ARGS "testinput/obserrorcovariance_test.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_hofx
                  TYPE SCRIPT
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_hofx.x testinput/hofx_makeobs3d.yml"             
                  DEPENDS soca_hofx.x )

ecbuild_add_test( TARGET test_soca_3dvar
                  TYPE SCRIPT
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_3dvar.x testinput/3dvar_test.yml"                       testoutput/3dvar.test
                  DEPENDS soca_3dvar.x )

ecbuild_add_test( TARGET test_soca_3dvarfgat
                  TYPE SCRIPT
                  COMMAND "./run.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_3dvar.x testinput/3dvarfgat_test.yml"
                  DEPENDS soca_3dvar.x )		  

ecbuild_add_test( TARGET test_soca_4denvar
                  TYPE SCRIPT
                  COMMAND "./compare.sh"
                  ARGS "mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/soca_3dvar.x testinput/4denvar_nicas_test.yml"
                  DEPENDS soca_3dvar.x )
		  
ecbuild_add_test( TARGET test_soca_dirac
                  TYPE SCRIPT
                  MPI     1		  
                  COMMAND "./compare.sh"
                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_dirac.x testinput/nicas_test.yml"
                       testoutput/dirac.test
                  DEPENDS soca_3dvar.x )		  
  
ecbuild_add_test( TARGET test_soca_balance
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/kst.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_bkgerr
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/bkgerr.yml"
                  LIBS    soca )

ecbuild_add_test( TARGET test_soca_vertconv
                  BOOST
                  SOURCES executables/TestVariableChange.cc
                  ARGS "testinput/vertconv.yml"
                  LIBS    soca )		
		
#ecbuild_add_test( TARGET test_soca_parameters_bump_cov
#                  TYPE SCRIPT
#                  COMMAND ${CMAKE_BINARY_DIR}/bin/compare.sh
#                  ARGS "${CMAKE_BINARY_DIR}/bin/soca_parameters.x testinput/parametersbump.cov.yml" )
