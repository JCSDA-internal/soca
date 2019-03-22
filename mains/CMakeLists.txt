ecbuild_add_executable( TARGET  soca_dirac.x
                        SOURCES socaDirac.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_forecast.x
                        SOURCES Forecast.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_makeobs.x
                        SOURCES MakeObs.cc
                        LIBS    soca
                      )


ecbuild_add_executable( TARGET  soca_hofx.x
                        SOURCES HofX.cc
                        LIBS    soca
                      )

ecbuild_add_executable( TARGET  soca_3dvar.x
                        SOURCES 3DVar.cc
                        LIBS    soca
                      )
		      
#ecbuild_add_executable( TARGET  soca_parameters.x
#                        SOURCES EstimateParams.cc
#                        LIBS    soca
#                      )

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

ecbuild_add_executable( TARGET  soca_checkpoint_model.x
                        SOURCES CheckpointModel.cc
                        LIBS    soca
                      )

# add soca coding norms test
ecbuild_add_test( TARGET soca_mains_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_BINARY_DIR}/bin/cpplint.py
                  ARGS --quiet --recursive ${CMAKE_CURRENT_SOURCE_DIR} )
