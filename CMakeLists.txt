cmake_minimum_required(VERSION 3.2)
enable_language(Fortran)

project(NS_DG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wextra -Wno-unused -pg")

set(SOURCE_FILES driver.cpp src/DG_Field_2d/DG_Field_2d.h src/DG_Field_2d/DG_Field_2d.cpp src/DG_Element_2d/DG_Element_2d.h src/DG_Element_2d/DG_Element_2d.cpp)

add_executable(nsdg ${SOURCE_FILES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})


find_package( BLAS REQUIRED )
message( STATUS BLAS found: ${BLAS_LIBRARIES} )
find_package( LAPACK REQUIRED )
message( STATUS LAPACK found: ${LAPACK_LIBRARIES} )

include_directories( ${BLAS_INCLUDE_DIRS}  ${LAPACK_INCLUDE_DIRS}  )

target_link_libraries(nsdg ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )