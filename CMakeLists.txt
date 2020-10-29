cmake_minimum_required(VERSION 3.0.0)
project(LDPC-with-SIMD VERSION 0.1.0)

include(CTest)
enable_testing()
set(SOURCES src/main.cpp src/alist/alist_matrix.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
