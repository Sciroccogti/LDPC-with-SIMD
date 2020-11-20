cmake_minimum_required(VERSION 3.0.0)
project(LDPC-with-SIMD VERSION 0.1.0)

# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11") set (CMAKE_CXX_FLAGS
# "${CMAKE_CXX_FLAGS} -std=c++11")

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
# if(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
# else()
#   check_cxx_compiler_flag("/arch:AVX" COMPILER_OPT_ARCH_AVX_SUPPORTED)
#   if(COMPILER_OPT_ARCH_AVX_SUPPORTED)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
#   endif()
# endif()

include(CTest)
enable_testing()

# Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

set(SOURCES src/main.cpp src/Alist/alist_matrix.c src/MatrixMath/MatrixMath.cpp
            src/LDPC/LDPC.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME}
                           PRIVATE ${PROJECT_SOURCE_DIR}/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
