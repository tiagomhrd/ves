cmake_minimum_required(VERSION 3.23.1)

project(tests)

add_executable(${PROJECT_NAME} 
${CMAKE_CURRENT_SOURCE_DIR}/src/Tests.cpp 
)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

target_include_directories(${PROJECT_NAME} 
PRIVATE src  ../ves/src ../third_party/Catch2/ ../third_party/ 
)

target_link_libraries(${PROJECT_NAME} 
Catch2 
ves 
)