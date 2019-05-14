# Config files for soca
list( APPEND soca_scripts
  cycle.sh
  cycle.py  
  )

list( APPEND soca_exp_input
  socainput/3dvar.yml
  socainput/3dvarfgat.yml
  socainput/4dhyb.yml
  socainput/cycle.yml
  socainput/genenspert.yml
  socainput/static_SocaError_init.yml
  socainput/static_SocaLoc_init.yml
)

list( APPEND soca_model_param
  ../test/Data/360x210x75/input-ymd.nml
  ../test/Data/360x210x75/MOM_input
  ../test/Data/360x210x75/MOM_override
  ../test/Data/360x210x75/data_table
  ../test/Data/360x210x75/diag_table
  ../test/Data/360x210x75/field_table
  ../test/Data/360x210x75/ice.bkgerror.nc
  ../test/Data/360x210x75/ocn.bkgerror.nc
  ../test/Data/rossrad.dat
  ../test/Data/godas_sst_bgerr.nc  
)

list( APPEND soca_model_restarts
  ../test/Data/360x210x75/INPUT/MOM.res.nc
  ../test/Data/360x210x75/INPUT/ice_model.res.nc
  ../test/Data/360x210x75/INPUT/layer_coord.nc
  ../test/Data/360x210x75/INPUT/hycom1_75_800m.nc
  ../test/Data/360x210x75/INPUT/ocean_hgrid.nc
  ../test/Data/360x210x75/INPUT/topog.nc
)


list( APPEND soca_model_forcing
  ../test/Data/360x210x75/INPUT/geothermal_heating_cm2g.nc
  ../test/Data/360x210x75/INPUT/seawifs_1998-2006_GOLD_smoothed_2X.nc
  ../test/Data/360x210x75/INPUT/tideamp.nc
  ../test/Data/360x210x75/INPUT/topog.nc
  ../test/Data/360x210x75/INPUT/gustiness_qscat.nc
  ../test/Data/360x210x75/INPUT/ocean_forcing_daily.nc
  ../test/Data/360x210x75/INPUT/ocean_precip_monthly.nc
  ../test/Data/360x210x75/INPUT/sgs_h2.nc
)

# Copy ioda converters
foreach(FILENAME ${soca_obs2ioda})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_BINARY_DIR}/bin/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

# Copy scripts
foreach(FILENAME ${soca_scripts})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

# Copy scripts
foreach(FILENAME ${soca_scripts})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

# Create directory for observations, DA output, ...
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

# Create directory for experiment input
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/expinput)
foreach(FILENAME ${soca_exp_input})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND}  -E create_symlink #shitch back to copy when debugging is done -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/expinput/${filename} )
endforeach(FILENAME)

# Create directory for model input and copy all files
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/INPUT)
foreach(FILENAME ${soca_model_restarts})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/INPUT/${filename} )
endforeach(FILENAME)

# MOM's resource files
foreach(FILENAME ${soca_model_param})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

# MOM's forcing files
foreach(FILENAME ${soca_model_forcing})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/INPUT/${filename} )
endforeach(FILENAME)

# Create MOM's output directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RESTART)

# Create directory for BUMP output
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bump)

# Test cycling

ecbuild_add_test( TARGET exp_soca_cycle
                  ENABLED OFF
                  TYPE SCRIPT
                  COMMAND "./cycle.py"
                  ARGS "./expinput/cycle.yml")
