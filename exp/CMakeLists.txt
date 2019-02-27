# Config files for soca
list( APPEND soca_scripts
  cycle.sh
  cycle.py  
  )

list( APPEND soca_exp_input
  socainput/3dvar.yml
  socainput/3dvarfgat.yml
  socainput/4dhyb.yml  
  socainput/genenspert.yml
  socainput/static_SocaError_init.yml
  socainput/static_SocaLoc_init.yml
)

list( APPEND soca_model_param
  ../test/Data/360x210x63/input-ymd.nml
  ../test/Data/360x210x63/MOM_input
  ../test/Data/360x210x63/MOM_override
  ../test/Data/360x210x63/data_table
  ../test/Data/360x210x63/diag_table
  ../test/Data/360x210x63/field_table
  ../test/Data/360x210x63/Vertical_coordinate.nc
  ../test/Data/360x210x63/ice.bkgerror.nc
  ../test/Data/360x210x63/ocn.bkgerror.nc
  ../test/Data/rossrad.dat
  ../test/Data/godas_sst_bgerr.nc  
)

list( APPEND soca_model_restarts
  ../test/Data/360x210x63/INPUT/MOM.res.nc
  ../test/Data/360x210x63/INPUT/ice_model.res.nc
  ../test/Data/360x210x63/INPUT/ocean_hgrid.nc
  ../test/Data/360x210x63/INPUT/topog.nc
)


list( APPEND soca_model_forcing
  ../test/Data/360x210x63/INPUT/geothermal_heating_cm2g.nc
  ../test/Data/360x210x63/INPUT/seawifs_1998-2006_GOLD_smoothed_2X.nc
  ../test/Data/360x210x63/INPUT/tideamp.nc
  ../test/Data/360x210x63/INPUT/topog.nc
  ../test/Data/360x210x63/INPUT/OM3_zgrid.nc
  ../test/Data/360x210x63/INPUT/gustiness_qscat.nc
  ../test/Data/360x210x63/INPUT/ocean_forcing_daily.nc
  ../test/Data/360x210x63/INPUT/ocean_precip_monthly.nc
  ../test/Data/360x210x63/INPUT/sgs_h2.nc
)

# Copy scripts
  message("================" ${soca_scripts})  
foreach(FILENAME ${soca_scripts})
  get_filename_component(filename ${FILENAME} NAME )
  message("================" $FILENAME)  
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)

# Create directory for observations, DA output, ...
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data)

# Create directory for experiment input
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/expinput)
foreach(FILENAME ${soca_exp_input})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
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
