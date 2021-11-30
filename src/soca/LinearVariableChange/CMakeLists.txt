add_subdirectory(Balance)
add_subdirectory(Base)
add_subdirectory(BkgErr)
add_subdirectory(BkgErrFilt)
add_subdirectory(BkgErrGodas)
add_subdirectory(HorizFilt)
add_subdirectory(Model2GeoVaLs)
add_subdirectory(util)
add_subdirectory(VertConv)

soca_target_sources(
  LinearVariableChange.cc
  LinearVariableChange.h
)
