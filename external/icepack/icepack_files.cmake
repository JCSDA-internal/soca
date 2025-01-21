# NOTE: only the files needed by SOCA have been enabled

#Icepack List:
list(APPEND icepack_files
#   Icepack/columnphysics/icepack_aerosol.F90
#   Icepack/columnphysics/icepack_age.F90
#   Icepack/columnphysics/icepack_algae.F90
#   Icepack/columnphysics/icepack_atmo.F90
#   Icepack/columnphysics/icepack_brine.F90
#   Icepack/columnphysics/icepack_firstyear.F90
#   Icepack/columnphysics/icepack_flux.F90
#   Icepack/columnphysics/icepack_fsd.F90
#   Icepack/columnphysics/icepack_intfc.F90
#   Icepack/columnphysics/icepack_isotope.F90
  Icepack/columnphysics/icepack_itd.F90
  Icepack/columnphysics/icepack_kinds.F90
#   Icepack/columnphysics/icepack_mechred.F90
#   Icepack/columnphysics/icepack_meltpond_cesm.F90
#   Icepack/columnphysics/icepack_meltpond_lvl.F90
#   Icepack/columnphysics/icepack_meltpond_topo.F90
  Icepack/columnphysics/icepack_mushy_physics.F90
#   Icepack/columnphysics/icepack_ocean.F90
#   Icepack/columnphysics/icepack_orbital.F90
  Icepack/columnphysics/icepack_parameters.F90
#   Icepack/columnphysics/icepack_shortwave.F90
#   Icepack/columnphysics/icepack_therm_0layer.F90
#   Icepack/columnphysics/icepack_therm_bl99.F90
#   Icepack/columnphysics/icepack_therm_itd.F90
#   Icepack/columnphysics/icepack_therm_mushy.F90
  Icepack/columnphysics/icepack_therm_shared.F90
#   Icepack/columnphysics/icepack_therm_vertical.F90
  Icepack/columnphysics/icepack_tracers.F90
  Icepack/columnphysics/icepack_warnings.F90
#   Icepack/columnphysics/icepack_wavefracspec.F90
#   Icepack/columnphysics/icepack_zbgc.F90
  Icepack/columnphysics/icepack_zbgc_shared.F90
#   Icepack/columnphysics/icepack_zsalinity.F90
)