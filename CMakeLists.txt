cmake_minimum_required(VERSION 2.8)

project(geolines)

set(MFEM_DIR /u/artemyev/software/mfem-3.0.1)
include_directories(${MFEM_DIR})
link_directories(${MFEM_DIR})

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} mfem)
