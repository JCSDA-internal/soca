list( APPEND soca_scripts
  marine2ioda.py
)

# Link scripts
foreach(FILENAME ${soca_scripts})
  get_filename_component(filename ${FILENAME} NAME )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${filename} )
endforeach(FILENAME)
