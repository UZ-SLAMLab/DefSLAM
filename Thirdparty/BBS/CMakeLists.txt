cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 14)
project(BBS VERSION 0.1 LANGUAGES CXX)
find_package(Eigen3 CONFIG REQUIRED)
IF (WIN32)
    set(LIB_TYPE)
ELSE()
    set(LIB_TYPE SHARED)
ENDIF()



add_library(${PROJECT_NAME} ${LIB_TYPE}
     bbs.cc
     bbs_coloc.cc
     )
target_include_directories(  ${PROJECT_NAME}
PUBLIC
    ${PROJECT_SOURCE_DIR} 
    
)

target_link_libraries(${PROJECT_NAME} ${EIGEN3_LIBS})
