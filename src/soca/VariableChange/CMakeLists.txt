add_subdirectory(Model2Ana)
add_subdirectory(Base)
add_subdirectory(Model2GeoVaLs)

if( ${icepack_FOUND} )
  add_subdirectory(Soca2Cice)
endif()

soca_target_sources(
  VariableChange.cc
  VariableChange.h
)
