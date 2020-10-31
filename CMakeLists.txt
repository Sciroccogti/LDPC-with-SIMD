cmake_minimum_required(VERSION 3.0.0)
project(LDPC-with-SIMD VERSION 0.1.0)

# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include(CTest)
enable_testing()

# Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

set(SOURCES src/main.cpp src/Alist/alist_matrix.c src/MatrixMath/MatrixMath.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
