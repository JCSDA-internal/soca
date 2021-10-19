
# Create symlink to MOM6 and soca static files
execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/../test/Data
           ${CMAKE_CURRENT_BINARY_DIR}/Data
           )

# Link tutorial scripts
set( tutorial_scripts
  tutorial_gridgen.sh
  tutorial_bump_op.sh
  tutorial_perturb_ic.sh
  tutorial_make_observations.sh
  tutorial_3dvar.sh
  tutorial_3dvarfgat.sh
  tutorial_synthetic_observations.py
  tutorial_plot.py
  tutorial_tools.sh
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
