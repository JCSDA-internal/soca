ecbuild_add_executable( TARGET  soca_forecast.x
                        SOURCES Forecast.cc
                        LIBS    ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS   ${UFS_LIBS}
)

ecbuild_add_executable( TARGET  soca_hofx.x
                        SOURCES HofX.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_enshofx.x
                        SOURCES EnsHofX.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_var.x
                        SOURCES Var.cc
                        LIBS   ${UFS_LIBS}
                                saber
                      )

ecbuild_add_executable( TARGET  soca_error_covariance_toolbox.x
                        SOURCES ErrorCovarianceToolbox.cc
                        LIBS   ${UFS_LIBS}
                                saber
                      )

ecbuild_add_executable( TARGET  soca_enspert.x
                        SOURCES GenEnsPertB.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_sqrtvertloc.x
                        SOURCES SqrtOfVertLoc.cc
                        LIBS    ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_ensmeanandvariance.x
                        SOURCES EnsMeanAndVariance.cc
                        LIBS    ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_ensrecenter.x
                        SOURCES EnsRecenter.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_hybridgain.x
                        SOURCES HybridGain.cc
                        LIBS   ${UFS_LIBS}
                      )

#ecbuild_add_executable( TARGET  soca_checkpoint_model.x
#                        SOURCES CheckpointModel.cc
#                        LIBS   ${UFS_LIBS}
#                      )

ecbuild_add_executable( TARGET  soca_gridgen.x
                        SOURCES GridGen.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_letkf.x
                        SOURCES LETKF.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_convertincrement.x
                        SOURCES ConvertIncrement.cc
                        LIBS    ${UFS_LIBS}
                        )

ecbuild_add_executable( TARGET  soca_convertstate.x
                        SOURCES ConvertState.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_addincrement.x
                        SOURCES AddIncrement.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_diffstates.x
                        SOURCES DiffStates.cc
                        LIBS   ${UFS_LIBS}
                      )

ecbuild_add_executable( TARGET  soca_setcorscales.x
                        SOURCES SetCorScales.cc
                        LIBS    ${UFS_LIBS}
                      )
