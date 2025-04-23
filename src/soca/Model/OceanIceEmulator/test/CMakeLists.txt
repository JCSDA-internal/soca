if(Torch_FOUND)
  ecbuild_add_executable( TARGET  soca_test_oceaniceffnn.x
                          SOURCES test_OceanIceFFNN.cc
                          LIBS    soca
                        )
endif()
