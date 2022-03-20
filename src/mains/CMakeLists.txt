ecbuild_add_executable( TARGET  soca_dirac.x
                        SOURCES Dirac.cc
                        LIBS    ${mylibs}
                      )

ecbuild_add_executable( TARGET  test_soca_model.x
                        SOURCES TestModel.cc
                        LIBS    ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_forecast.x
                        SOURCES Forecast.cc
                        LIBS    
                         soca
                         cdeps
                         cmeps
                         piof
                         pioc
                         stochastic_physics
                         cice
                         esmf
                      )

ecbuild_add_executable( TARGET  soca_hofx3d.x
                        SOURCES HofX3D.cc
                        LIBS  ${mylibs}
)

ecbuild_add_executable( TARGET  soca_hofx.x
                        SOURCES HofX.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_enshofx.x
                        SOURCES EnsHofX.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_var.x
                        SOURCES Var.cc
                        LIBS   ${mylibs}
                                saber
                      )

ecbuild_add_executable( TARGET  soca_error_covariance_training.x
                        SOURCES ErrorCovarianceTraining.cc
                        LIBS   ${mylibs}
                                saber
                      )

ecbuild_add_executable( TARGET  soca_staticbinit.x
                        SOURCES StaticBInit.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_enspert.x
                        SOURCES GenEnsPertB.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_sqrtvertloc.x
                        SOURCES SqrtOfVertLoc.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_ensvariance.x
                        SOURCES EnsVariance.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_ensrecenter.x
                        SOURCES EnsRecenter.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_hybridgain.x
                        SOURCES HybridGain.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_checkpoint_model.x
                        SOURCES CheckpointModel.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_gridgen.x
                        SOURCES GridGen.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_letkf.x
                        SOURCES LETKF.cc
                        LIBS   ${mylibs}
                      )
ecbuild_add_executable( TARGET  soca_convertstate.x
                        SOURCES ConvertState.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_addincrement.x
                        SOURCES AddIncrement.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_diffstates.x
                        SOURCES DiffStates.cc
                        LIBS   ${mylibs}
                      )

ecbuild_add_executable( TARGET  soca_setcorscales.x
                        SOURCES SetCorScales.cc
                        LIBS    ${mylibs}
                      )
