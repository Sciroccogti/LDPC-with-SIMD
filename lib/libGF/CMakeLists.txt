cmake_minimum_required(VERSION 3.12.0)
project(libGF)

set(SOURCES src/GF.cpp)

set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CXX)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

target_include_directories(${PROJECT_NAME}
                           PRIVATE ${PROJECT_SOURCE_DIR}/include)
