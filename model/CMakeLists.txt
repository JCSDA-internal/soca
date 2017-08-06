list( APPEND oops_mom5cice5_src_files
Fields.cc
Fields.h
Fortran.h
Geometry.cc
Geometry.h
Increment.cc
Increment.h
State.cc
State.h
Traits.h
Variables.cc
Variables.h
mom5cice5_constants.f90
mom5cice5_fields.F90
mom5cice5_fields_interface.F90
mom5cice5_geom_mod.F90
mom5cice5_invent_state.f90
mom5cice5_vars_mod.F90
mom5cice5_ncutils_mod.F90
)

if ( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  set( cfdealloc C_F_NO_DEALLOC )
endif()

ecbuild_add_library( TARGET   mom5cice5
                     DEFINITIONS ${cfdealloc}
                     SOURCES  ${oops_mom5cice5_src_files}
                     LIBS     oops ${LAPACK_LIBRARIES}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )

