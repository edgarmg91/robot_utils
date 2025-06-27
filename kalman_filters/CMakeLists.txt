cmake_minimum_required(VERSION 2.8.3)
project(kalman_filter)

#Compile witch C++11
add_compile_options(-std=c++11)


#Setup Eigen
set(Eigen_DIR "/usr/include/eigen3")

## Node headers
include_directories(
  include 
  ${Eigen_DIR}
)

#Human-Pose node
add_executable(
  kalman_cpp
  src/main.cpp
  include/kalman_filter.cpp
)

target_link_libraries(
  kalman_cpp
)
