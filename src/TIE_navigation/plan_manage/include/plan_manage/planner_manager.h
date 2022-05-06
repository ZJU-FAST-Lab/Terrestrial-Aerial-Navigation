/**
* This file is part of Fast-Planner.
*
* Copyright 2019 Boyu Zhou, Aerial Robotics Group, Hong Kong University of Science and Technology, <uav.ust.hk>
* Developed by Boyu Zhou <bzhouai at connect dot ust dot hk>, <uv dot boyuzhou at gmail dot com>
* for more information see <https://github.com/HKUST-Aerial-Robotics/Fast-Planner>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* Fast-Planner is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Fast-Planner is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Fast-Planner. If not, see <http://www.gnu.org/licenses/>.
*/



#ifndef _PLANNER_MANAGER_H_
#define _PLANNER_MANAGER_H_

#include <bspline_opt/bspline_optimizer.h>
#include <bspline/non_uniform_bspline.h>
#include <uav_utils/geometry_utils.h>
#include <path_searching/astar.h>
#include <path_searching/kinodynamic_astar.h>
#include <plan_env/edt_environment.h>

#include <plan_manage/plan_container.hpp>

#include <ros/ros.h>

#include <traj_utils/polynomial_traj.h>

namespace fast_planner {

// Fast Planner Manager
// Key algorithms of mapping and planning are called

class FastPlannerManager {
  // SECTION stable
public:
  FastPlannerManager();
  ~FastPlannerManager();

  /* main planning interface */
  bool kinodynamicReplan(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel, Eigen::Vector3d start_acc,
                         Eigen::Vector3d end_pt, Eigen::Vector3d end_vel);
                         
  void initPlanModules(ros::NodeHandle& nh);
  void setGlobalWaypoints(vector<Eigen::Vector3d>& waypoints);

  int checkTrajCollision(double& distance, double& duration);
  inline int getPlanningType(){
    return planning_type_;
  }
  PlanParameters pp_;
  LocalTrajData local_data_;
  MidPlanData plan_data_;
  EDTEnvironment::Ptr edt_environment_;
  PolynomialTraj best_traj;
  double ground_judge;

private:
  /* main planning algorithms & modules */
  SDFMap::Ptr sdf_map_;
  bool close_goal_traj_;
  unique_ptr<Astar> geo_path_finder_;
  unique_ptr<KinodynamicAstar> kino_path_finder_;
  vector<BsplineOptimizer::Ptr> bspline_optimizers_;
  int continous_failures_count_{0};
  int primitive_num_;
  double search_time, optimization_time;
  int plan_round;
  int planning_type_;
  double odom_yaw;
  double replan_thresh_;
  double geo_astar_resolution_;
  
  Eigen::Vector3d odom_pos_, odom_vel_;  // odometry state
  Eigen::Quaterniond odom_orient_;

  ros::Subscriber odom_sub_;
  void odometryCallback(const nav_msgs::OdometryConstPtr& msg);
  void updateTrajInfo();

  // !SECTION stable
  // SECTION developing

public:
  typedef unique_ptr<FastPlannerManager> Ptr;

  // !SECTION
};
}  // namespace fast_planner

#endif