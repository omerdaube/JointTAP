cmake_minimum_required(VERSION 3.21)
project(JointTAP)

set(CMAKE_CXX_STANDARD 23)
include_directories(${CMAKE_SOURCE_DIR}/include)
set(CMAKE_CXX_FLAGS "-O0")

add_executable(JointTAP main.cpp utils.cpp incremental_otp.cpp joint_tap.cpp lp_bound.cpp assistance_otp.cpp)