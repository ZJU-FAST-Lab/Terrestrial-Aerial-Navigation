cmake_minimum_required(VERSION 2.8.3)
project(traj_utils)

find_package(catkin REQUIRED COMPONENTS
  bspline
  roscpp
  std_msgs
  cv_bridge
)

find_package(Eigen3 REQUIRED)
find_package(PCL 1.7 REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})

find_path(EIGEN_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen3
	/usr/include/eigen3
)

catkin_package(
 INCLUDE_DIRS include
 LIBRARIES traj_utils
 CATKIN_DEPENDS bspline
#  DEPENDS system_lib
)

include_directories( 
    SYSTEM 
    include 
    ${catkin_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
)

link_directories(${PCL_LIBRARY_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -O3 -Wall")

add_library( traj_utils
    src/planning_visualization.cpp 
    src/poly_traj_utils.cpp 
    src/polynomial_traj.cpp
    )
target_link_libraries( traj_utils
    ${catkin_LIBRARIES} 
    )  
