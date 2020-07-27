cmake_minimum_required(VERSION 2.8)
project(ToolsPCL)

find_package(PCL 1.5 REQUIRED COMPONENTS io common features surface)
message(STATUS "Using PCL ${PCL_VERSION}")
include_directories(${PROJECT_SOURCE_DIR}
			${PCL_INCLUDE_DIRS})

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_library(${PROJECT_NAME} SHARED
    PCLNormalEstimator.cc
    SmootherMLS.cc
)
target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES})


set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
