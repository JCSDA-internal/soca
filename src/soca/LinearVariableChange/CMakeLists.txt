add_subdirectory(Balance)
add_subdirectory(Base)
add_subdirectory(BkgErr)
add_subdirectory(BkgErrFilt)
add_subdirectory(BkgErrGodas)
add_subdirectory(LinearModel2GeoVaLs)
add_subdirectory(util)

soca_target_sources(
  LinearVariableChange.cc
  LinearVariableChange.h
)
