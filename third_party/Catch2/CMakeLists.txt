cmake_minimum_required(VERSION 3.23.1)

project(Catch2)

add_library(${PROJECT_NAME} 
${CMAKE_CURRENT_SOURCE_DIR}/catch_amalgamated.cpp
)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)