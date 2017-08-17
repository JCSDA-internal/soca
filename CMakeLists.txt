cmake_minimum_required(VERSION 3.6)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/model )

add_subdirectory( mains )
add_subdirectory( model )
add_subdirectory( test )

