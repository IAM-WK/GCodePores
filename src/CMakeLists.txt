cmake_minimum_required(VERSION 2.8)
 
project(GCodePores)

#SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "C:/Libs/boost_1_70_0")
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:/Libs/boost_1_70_0/stage/lib")

FIND_PACKAGE(Boost 1.70.0 COMPONENTS program_options)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

add_executable(GCodePores GCodePores.cpp helperfunctions.cpp helperfunctions.h GCodeParser.cpp GCodeParser.h GCodeTransform.cpp GCodeTransform.h GCodeInterpolation.cpp GCodeInterpolation.h PathBase.cpp PathBase.h mhdParser.cpp mhdParser.h GCodeAnalysis.cpp GCodeAnalysis.h TimeMeasurement.cpp TimeMeasurement.h PoreParser.cpp PoreParser.h )

