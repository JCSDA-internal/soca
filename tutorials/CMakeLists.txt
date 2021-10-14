
# Create symlink to MOM6 and soca static files
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/../test/Data
           ${CMAKE_CURRENT_BINARY_DIR}/Data
           )

# Link tutorial scripts
set( tutorial_scripts
  input.nml.sh
  3dvarfgat.yaml.sh
  checkpoint.yaml.sh
  forecast.yaml.sh
  prep.mom6-soca.static.sh
  soca_3dvarfgat.sh
  soca_ensemble_perturbations.sh
  soca_gridgen.sh
  soca_plot.py
  soca_static_correlation_op.sh
  soca_synthetic_observations.py
  soca_synthetic_observations.sh
)
foreach(FILENAME ${tutorial_scripts})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

# Create symlink to tutorial yaml configuration files
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/config
           ${CMAKE_CURRENT_BINARY_DIR}/config
           )

# Create symlink to JEDI executables
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_BINARY_DIR}/../../bin
           ${CMAKE_CURRENT_BINARY_DIR}/bin
           )
