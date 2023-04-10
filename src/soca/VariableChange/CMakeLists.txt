add_subdirectory(Ana2Model)
add_subdirectory(Base)
add_subdirectory(Model2GeoVaLs)

if ( BUILD_ICEPACK )
  add_subdirectory(Soca2Cice)
endif()

soca_target_sources(
  VariableChange.cc
  VariableChange.h
)
