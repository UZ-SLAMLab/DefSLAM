cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 14)
project(ORBSLAM VERSION 0.1 LANGUAGES CXX)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(OpenCV 4.0.0 REQUIRED)
find_package(Pangolin REQUIRED)
FIND_PACKAGE(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)

IF (WIN32)
    add_compile_options("/openmp:experimental")
    set(LIB_TYPE)
ELSE()
    set(LIB_TYPE SHARED)
ENDIF()

include_directories(
${PROJECT_SOURCE_DIR}/ThirdParty/ORBSLAM2
${EIGEN3_INCLUDE_DIR}
${Pangolin_INCLUDE_DIRS}
)

FILE(GLOB SRC_ORBSLAM2_FILES
     "src/*.cc") 

add_library(${PROJECT_NAME}  ${LIB_TYPE}
${SRC_ORBSLAM2_FILES}
)
IF (WIN32)
target_compile_options(${PROJECT_NAME} PRIVATE /bigobj)
endif()
target_link_libraries(${PROJECT_NAME}
${OpenCV_LIBS}
${EIGEN3_LIBS} 
${Pangolin_LIBRARIES}
${CERES_LIBRARIES}
DBoW2
g2o
)

