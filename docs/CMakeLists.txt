# first we can indicate the documentation build as an option and
# set it to OFF by default since documenation will usually be
# built at the bundle level
option(ENABLE_SOCA_DOC "Build SOCA documentation" OFF)

if (ENABLE_SOCA_DOC)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        message(STATUS "Build SOCA Documentation")
        include(GNUInstallDirs)

        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

        add_custom_target( doc_doxygen_soca ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating SOCA documentation with Doxygen"
            VERBATIM )
        set( HTML_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/html )
        install(DIRECTORY ${HTML_BUILD_DIR} DESTINATION ${CMAKE_INSTALL_FULL_DATADIR}/${PROJECT_NAME}/doc COMPONENT Documentation)
    else ()
        message(WARNING "Doxygen need to be installed to generate SOCA documentation")
    endif ()
endif ()
