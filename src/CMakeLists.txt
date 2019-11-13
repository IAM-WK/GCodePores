cmake_minimum_required(VERSION 2.8)
 
project(GCodePores)

# try to include boost
#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "C:/Libs/boost_1_67_0")
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:/Libs/boost_1_67_0/stage/lib")

#FIND_PACKAGE(Boost 1.67.0 COMPONENTS program_options)
#IF (Boost_FOUND)
#    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
#    ADD_DEFINITIONS( "-DHAS_BOOST" )
#ENDIF()
 # end try to include boost
 
#find_package(ITK REQUIRED)
#include(${ITK_USE_FILE})
#if (ITKVtkGlue_LOADED)
#  find_package(VTK REQUIRED)
#  include(${VTK_USE_FILE})
#endif()
 
add_executable(GCodePores GCodePores.cpp helperfunctions.cpp helperfunctions.h GCodeParser.cpp GCodeParser.h GCodeTransform.cpp GCodeTransform.h GCodeInterpolation.cpp GCodeInterpolation.h PathBase.cpp PathBase.h)
 
#if( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(GCodeImageDirection ITKReview ${ITK_LIBRARIES})
#else( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(GCodeImageDirection ${ITK_LIBRARIES})
#endif( "${ITK_VERSION_MAJOR}" LESS 4 )