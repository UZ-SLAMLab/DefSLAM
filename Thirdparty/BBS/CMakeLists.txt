cmake_minimum_required(VERSION 2.8)
project(BBS)

find_package(Eigen3 REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR} 
    ${EIGEN3_INCLUDE_DIR}
)

add_library(${PROJECT_NAME} SHARED
     bbs.cc
     bbs_coloc.cc
     )

target_link_libraries(${PROJECT_NAME} ${EIGEN3_LIBS})
