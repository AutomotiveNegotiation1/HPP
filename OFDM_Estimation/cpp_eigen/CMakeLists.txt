cmake_minimum_required(VERSION 2.8)
project(OFDM_Estimation)

include_directories(
${PROJECT_SOURCE_DIR}
${PROJECT_SOURCE_DIR}/include
${PROJECT_SOURCE_DIR}/lib
)

# find_package(Boost REQUIRED) 
find_package(Eigen3 REQUIRED) 

if (Eigen3_FOUND)
    include_directories(${Eigen3_INCLUDE_DIRS})
    
else()
    message(FATAL_ERROR "Eigen3 not found")
endif()

#include_directories(${Boost_INCLUDE_DIRS})

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(output ${SOURCES})

#target_link_libraries(output ${Boost_LIBRARIES} ${Eigen3_LIBRARIES})
target_link_libraries(output ${Eigen3_LIBRARIES})