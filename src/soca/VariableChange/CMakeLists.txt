add_subdirectory(Base)
add_subdirectory(Model2Ana)
add_subdirectory(Model2GeoVaLs)
add_subdirectory(Soca2Cice)

soca_target_sources(
  VariableChange.cc
  VariableChange.h
)
