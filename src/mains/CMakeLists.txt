ecbuild_add_executable( TARGET  soca_forecast.x
                        SOURCES Forecast.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS    soca
)

ecbuild_add_executable( TARGET  soca_hofx.x
                        SOURCES HofX.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_enshofx.x
                        SOURCES EnsHofX.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_var.x
                        SOURCES Var.cc
                        LIBS    soca
                                saber
                      )

ecbuild_add_executable( TARGET  soca_error_covariance_toolbox.x
                        SOURCES ErrorCovarianceToolbox.cc
                        LIBS    soca
                                saber
                      )

ecbuild_add_executable( TARGET  soca_enspert.x
                        SOURCES GenEnsPertB.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_sqrtvertloc.x
                        SOURCES SqrtOfVertLoc.cc
                        LIBS    soca
                        )

ecbuild_add_executable( TARGET  soca_ensmeanandvariance.x
                        SOURCES EnsMeanAndVariance.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_ensrecenter.x
                        SOURCES EnsRecenter.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_hybridgain.x
                        SOURCES HybridGain.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_gridgen.x
                        SOURCES GridGen.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_letkf.x
                        SOURCES LETKF.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_convertincrement.x
                        SOURCES ConvertIncrement.cc
                        LIBS    soca
                        )

ecbuild_add_executable( TARGET  soca_convertstate.x
                        SOURCES ConvertState.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_addincrement.x
                        SOURCES AddIncrement.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_diffstates.x
                        SOURCES DiffStates.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_setcorscales.x
                        SOURCES SetCorScales.cc
                        LIBS    soca
                      )
ecbuild_add_executable( TARGET  soca_gen_hybrid_linear_model_coeffs.x
                        SOURCES GenHybridLinearModelCoeffs.cc
                        LIBS    soca
                      )
ecbuild_add_executable( TARGET  soca_linearization_error.x
                        SOURCES LinearizationError.cc
                        LIBS    soca
                      )

if(Torch_FOUND)
  ecbuild_add_executable( TARGET  soca_ice_emulator.x
                        SOURCES IceEmul.cc
                        LIBS    soca NetCDF::NetCDF_C
                      )
endif()