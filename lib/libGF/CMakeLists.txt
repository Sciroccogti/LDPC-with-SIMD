cmake_minimum_required(VERSION 3.12.0)
project(libGF)

set(SOURCES src/GF.cpp)

set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CXX)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

target_include_directories(${PROJECT_NAME}
                           PRIVATE ${PROJECT_SOURCE_DIR}/include)

# # GPROF
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

# auto optimize
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
