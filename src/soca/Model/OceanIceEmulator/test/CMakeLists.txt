if(Torch_FOUND)
  ecbuild_add_executable( TARGET  test_oceaniceffnn.x
                          SOURCES test_OceanIceFFNN.cc
                          LIBS    soca
                        )
endif()
