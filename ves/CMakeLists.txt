cmake_minimum_required(VERSION 3.23.1)

project(ves)

add_library(${PROJECT_NAME} STATIC
${CMAKE_CURRENT_SOURCE_DIR}/src/ves.h 
${CMAKE_CURRENT_SOURCE_DIR}/src/ves_internal.h 
${CMAKE_CURRENT_SOURCE_DIR}/src/ves_internal.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/src/V2D.h 
${CMAKE_CURRENT_SOURCE_DIR}/src/V2D.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/src/SV2D.h 
${CMAKE_CURRENT_SOURCE_DIR}/src/SV2D.cpp 
${CMAKE_CURRENT_SOURCE_DIR}/src/SE2V2D.h 
${CMAKE_CURRENT_SOURCE_DIR}/src/SE2V2D.cpp 
)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

target_include_directories(${PROJECT_NAME} 
PRIVATE src ../third_party/ 
)

target_link_libraries(${PROJECT_NAME} 
ptp 
)