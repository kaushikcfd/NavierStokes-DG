cmake_minimum_required(VERSION 3.2)

project(NS_DG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -llapacke -lblas")
SET(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wextra -Wno-unused -pg")

set(SOURCE_FILES driver.cpp driver.h src/DG_Field_2d/DG_Field_2d.h src/DG_Field_2d/DG_Field_2d.cpp src/DG_Element_2d/DG_Element_2d.h src/DG_Element_2d/DG_Element_2d.cpp)

add_executable(nsdg ${SOURCE_FILES})
