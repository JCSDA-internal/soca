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

#ecbuild_add_executable( TARGET  soca_moments.x
#                        SOURCES Moments.cc
#                        LIBS    soca
#                      )
