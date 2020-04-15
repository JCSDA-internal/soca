ecbuild_add_executable( TARGET  soca_dirac.x
                        SOURCES Dirac.cc
                        LIBS    soca
                                saber
                      )

ecbuild_add_executable( TARGET  soca_forecast.x
                        SOURCES Forecast.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_makeobs.x
                        SOURCES MakeObs.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_hofx_nomodel.x
                        SOURCES HofXNoModel.cc
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

ecbuild_add_executable( TARGET  soca_3dvar.x
                        SOURCES 3DVar.cc
                        LIBS    soca
                                saber
                      )

ecbuild_add_executable( TARGET  soca_parameters.x
                        SOURCES EstimateParams.cc
                        LIBS    soca
                                saber
                      )

ecbuild_add_executable( TARGET  soca_staticbinit.x
                        SOURCES StaticBInit.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_enspert.x
                        SOURCES GenEnsPertB.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_ensvariance.x
                        SOURCES EnsVariance.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_ensrecenter.x
                        SOURCES EnsRecenter.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_checkpoint_model.x
                        SOURCES CheckpointModel.cc
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
ecbuild_add_executable( TARGET  soca_convertstate.x
                        SOURCES ConvertState.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_addincrement.x
                        SOURCES AddIncrement.cc
                        LIBS    soca
                      )
