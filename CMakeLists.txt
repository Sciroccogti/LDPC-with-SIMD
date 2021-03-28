cmake_minimum_required(VERSION 3.12.0)
project(LDPC-with-SIMD VERSION 0.1.0)
include(CTest)
enable_testing()

add_subdirectory(${PROJECT_SOURCE_DIR}/lib/libGF)

set(CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11") set (CMAKE_CXX_FLAGS
# "${CMAKE_CXX_FLAGS} -std=c++11")

# # OPENMP find_package(OpenMP REQUIRED) if(OPENMP_FOUND) set(CMAKE_C_FLAGS
# "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}") set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}
# ${OpenMP_CXX_FLAGS}") set(CMAKE_EXE_LINKER_FLAGS
# "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

# pthread conflict with OPENMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

# SIMD
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
if(COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mavx -mfma")
  # else() check_cxx_compiler_flag("/arch:AVX" COMPILER_OPT_ARCH_AVX_SUPPORTED)
  # if(COMPILER_OPT_ARCH_AVX_SUPPORTED) set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}
  # /arch:AVX") endif()
endif()

# MIPP add_definitions("-DMIPP_ALIGNED_LOADS")
# add_definitions("-DMIPP_NO_INTRINSICS")

set(SOURCES
    src/main.cpp
    src/utils/utils.cpp
    src/Alist/Alist.cpp
    src/Alist/alist_matrix.c
    src/Alist/nbalist_matrix.c
    src/MatrixMath/MatrixMath.cpp
    src/MatrixMath/NBMatrixMath.cpp
    src/LDPC/LDPC.cpp
    src/LDPC/NBLDPC.cpp
    src/LDPC/Tanner.cpp
    src/Modem/Modem.cpp
    src/Channel/Channel.cpp
    src/UI/pyplot.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

# Eigen3
find_package(Eigen3 REQUIRED)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen3)

# YAML
find_package(yaml-cpp REQUIRED)
target_link_libraries(${PROJECT_NAME} yaml-cpp::yaml-cpp)

# Python
find_package(PythonLibs 3.6 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})

# GF
target_link_libraries(${PROJECT_NAME} libGF)
target_include_directories(${PROJECT_NAME}
                           PRIVATE ${PROJECT_SOURCE_DIR}/lib/libGF/include)

target_include_directories(${PROJECT_NAME}
                           PRIVATE ${PROJECT_SOURCE_DIR}/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
