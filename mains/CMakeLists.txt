ecbuild_add_executable( TARGET  mom5cice5_dirac.x
                        SOURCES mom5cice5Dirac.cc
                        LIBS    mom5cice5
                      )

ecbuild_add_executable( TARGET  mom5cice5_forecast.x
                        SOURCES Forecast.cc
                        LIBS    mom5cice5
                      )

