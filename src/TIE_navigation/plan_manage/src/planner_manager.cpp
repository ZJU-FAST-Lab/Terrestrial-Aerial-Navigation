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



// #include <fstream>
#include <plan_manage/planner_manager.h>
#include <thread>

namespace fast_planner {

// SECTION interfaces for setup and query

FastPlannerManager::FastPlannerManager() {}

FastPlannerManager::~FastPlannerManager() { std::cout << "des manager" << std::endl; }

void FastPlannerManager::initPlanModules(ros::NodeHandle& nh) {
  /* read algorithm parameters */

  nh.param("manager/max_vel", pp_.max_vel_, -1.0);
  nh.param("manager/max_acc", pp_.max_acc_, -1.0);
  nh.param("manager/max_jerk", pp_.max_jerk_, -1.0);
  nh.param("manager/dynamic_environment", pp_.dynamic_, -1);
  nh.param("manager/clearance_threshold", pp_.clearance_, -1.0);
  nh.param("manager/local_segment_length", pp_.local_traj_len_, -1.0);
  nh.param("manager/control_points_distance", pp_.ctrl_pt_dist, -1.0);

  bool use_geometric_path, use_kinodynamic_path, use_optimization, use_active_perception;
  nh.param("manager/use_geometric_path", use_geometric_path, true);
  nh.param("manager/use_kinodynamic_path", use_kinodynamic_path, true);
  nh.param("manager/use_optimization", use_optimization, true);
  nh.param("manager/planning_type", planning_type_, 0);
  nh.param("sdf_map/ground_judge", ground_judge, 1.0);
  nh.param("manager/primitive_num", primitive_num_, 1);
  nh.param("fsm/thresh_replan", replan_thresh_, -1.0);
  nh.param("astar/resolution_astar", geo_astar_resolution_, -1.0);
  local_data_.traj_id_ = 0;
  local_data_.best_primitive_index_ = 0;
  sdf_map_.reset(new SDFMap);
  sdf_map_->initMap(nh);
  edt_environment_.reset(new EDTEnvironment);
  edt_environment_->setMap(sdf_map_);
  
  odom_sub_ = nh.subscribe("/odom_world", 1, &FastPlannerManager::odometryCallback, this);
  odom_yaw = 0;

  if (use_geometric_path) {
    geo_path_finder_.reset(new Astar);
    geo_path_finder_->setParam(nh);
    geo_path_finder_->setEnvironment(edt_environment_);
    geo_path_finder_->init();
    local_data_.motion_primitive_wpts_.resize(3 * primitive_num_);
  }

  if (use_kinodynamic_path) {
    kino_path_finder_.reset(new KinodynamicAstar);
    kino_path_finder_->setParam(nh);
    kino_path_finder_->setEnvironment(edt_environment_);
    kino_path_finder_->init();
  }

  if (use_optimization) {
    bspline_optimizers_.resize(10);
    for (int i = 0; i < 10; ++i) {
      bspline_optimizers_[i].reset(new BsplineOptimizer);
      bspline_optimizers_[i]->setParam(nh);
      bspline_optimizers_[i]->setEnvironment(edt_environment_);
    }
  }

  search_time = 0;
  optimization_time = 0;
  plan_round = 0;
}

void FastPlannerManager::setGlobalWaypoints(vector<Eigen::Vector3d>& waypoints) {
  plan_data_.global_waypoints_ = waypoints;
}
void FastPlannerManager::odometryCallback(const nav_msgs::OdometryConstPtr& msg) {
  odom_pos_(0) = msg->pose.pose.position.x;
  odom_pos_(1) = msg->pose.pose.position.y;
  odom_pos_(2) = msg->pose.pose.position.z;

  odom_vel_(0) = msg->twist.twist.linear.x;
  odom_vel_(1) = msg->twist.twist.linear.y;
  odom_vel_(2) = msg->twist.twist.linear.z;

  odom_orient_.w() = msg->pose.pose.orientation.w;
  odom_orient_.x() = msg->pose.pose.orientation.x;
  odom_orient_.y() = msg->pose.pose.orientation.y;
  odom_orient_.z() = msg->pose.pose.orientation.z;

  Eigen::Vector3d ypr = uav_utils::quaternion_to_ypr(odom_orient_);
  odom_yaw = ypr(0);
  // cout << "odom_yaw: " << odom_yaw <<endl;

}
int FastPlannerManager::checkTrajCollision(double& distance, double& duration) {

  int result = 0;
  double t_now = (ros::Time::now() - local_data_.start_time_).toSec();

  if(planning_type_ == 0){
    double tm, tmp;
    local_data_.position_traj_.getTimeSpan(tm, tmp);
    Eigen::Vector3d cur_pt = local_data_.position_traj_.evaluateDeBoor(tm + t_now);

    double          radius = 0.0;
    Eigen::Vector3d fut_pt;
    double          fut_t = 0.02;

    while (radius < 6.0 && t_now + fut_t < duration) {
      fut_pt = local_data_.position_traj_.evaluateDeBoor(tm + t_now + fut_t);

      double dist = edt_environment_->evaluateCoarseEDT(fut_pt, -1.0);
      if (dist < 0.1) {
        if(fut_t == 0){
          cout << "current quadrotor in collision!" << endl;
          result = 2;
          return result;
        }
        distance = radius;
        result = 1;
        return result;
      }

      radius = (fut_pt - cur_pt).norm();
      fut_t += 0.02;
    }
  }
  else{
    double duration = best_traj.getTimeSum();
    t_now = min(duration, t_now);
    Eigen::Vector3d cur_pt = best_traj.evaluate(t_now);

    double          radius = 0.0;
    Eigen::Vector3d fut_pt;
    double          fut_t = 0.02;
    while (radius < 6.0 && t_now + fut_t < duration - 1e-2) {
      fut_pt = best_traj.evaluate(t_now + fut_t);

      double dist = edt_environment_->evaluateCoarseEDT(fut_pt, -1.0);
      if (dist < 0.1) {
        if(fut_t == 0){
          cout << "current quadrotor in collision!" << endl;
          result = 2;
          return result;
        }
        result = 1;
        distance = radius;
        return result;
      }


      radius = (fut_pt - cur_pt).norm();
      fut_t += 0.02;
    }
  }


  return 0;
}

// !SECTION

bool FastPlannerManager::kinodynamicReplan(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel,
                                           Eigen::Vector3d start_acc, Eigen::Vector3d end_pt,
                                           Eigen::Vector3d end_vel) {

  std::cout << "[kino replan]: -----------------------" << std::endl;
  // cout << "start: " << start_pt.transpose() << ", " << start_vel.transpose() << ", "
  //      << start_acc.transpose() << "\ngoal:" << end_pt.transpose() << ", " << end_vel.transpose()
  //      << endl;

  if ((start_pt - end_pt).norm() < 0.2) {
    cout << "Close goal" << endl;
    continous_failures_count_++;
    return false;
  }

  ros::Time t1, t2;

  local_data_.start_time_ = ros::Time::now();
  double t_search = 0.0, t_opt = 0.0, t_adjust = 0.0;

  Eigen::Vector3d init_pos = start_pt;
  Eigen::Vector3d init_vel = start_vel;
  Eigen::Vector3d init_acc = start_acc;

  Eigen::Vector3d end_acc;
  end_acc.setZero();
  // kinodynamic path searching

  Eigen::MatrixXd ctrl_pts;
  double                  ts = pp_.ctrl_pt_dist / pp_.max_vel_;


  if(planning_type_ == 0){
    // proposed
    t1 = ros::Time::now();
    
    kino_path_finder_->reset();

    int status = kino_path_finder_->search(start_pt, start_vel, start_acc, end_pt, end_vel, false);

    if (status == KinodynamicAstar::NO_PATH) {
      cout << "[kino replan]: kinodynamic search fail!" << endl;

      // retry searching with discontinuous initial state
      kino_path_finder_->reset();
      status = kino_path_finder_->search(start_pt, start_vel, start_acc, end_pt, end_vel, false);

      if (status == KinodynamicAstar::NO_PATH) {
        cout << "[kino replan]: Can't find path." << endl;
        return false;
      } else {
        cout << "[kino replan]: retry search success." << endl;
      }

    } else {
      cout << "[kino replan]: kinodynamic search success." << endl;
    }

    t_search = (ros::Time::now() - t1).toSec();
    //for visualization
    plan_data_.kino_path_ = kino_path_finder_->getKinoTraj(0.01);
    vector<PathNodePtr> hybrid_path = kino_path_finder_->getGridPath();

    // parameterize the path to bspline
   
    vector<Eigen::Vector3d> hybrid_point_set, hybrid_start_end_derivatives, point_set, start_end_derivatives;
    vector<vector<Eigen::Vector3d>> hybrid_point_set_list; 
    vector<int> hybrid_motion_list, each_motion_list;

    kino_path_finder_->getHybridSamples(ts, hybrid_point_set_list, hybrid_start_end_derivatives, hybrid_motion_list);
    
    vector<vector<Eigen::Vector3d>> sampled_wpts;

    for(int i = 0; i < hybrid_point_set_list.size(); i++){
      for(int j = 0; j < hybrid_point_set_list[i].size(); j++){
        each_motion_list.push_back(hybrid_motion_list[i]);
        Eigen::Vector3d sampled_hybrid_pos = hybrid_point_set_list[i][j];
        hybrid_point_set.push_back(sampled_hybrid_pos);
      }
    }

    NonUniformBspline::parameterizeToBspline(ts, hybrid_point_set, hybrid_start_end_derivatives, ctrl_pts);
    NonUniformBspline init(ctrl_pts, 3, ts);
    // bspline trajectory optimization

    t1 = ros::Time::now();

    int cost_function = BsplineOptimizer::NORMAL_PHASE;

    if (status != KinodynamicAstar::REACH_END) {
      cost_function |= BsplineOptimizer::ENDPOINT;
    }
    if (hybrid_motion_list.size() > 1 || hybrid_motion_list[0] == 0){ 
      //Rolling/Ground state exists
      cost_function |= BsplineOptimizer::NON_HOLONOMIC;
    }
    bspline_optimizers_[0]->setMotionState(each_motion_list);
    ctrl_pts = bspline_optimizers_[0]->BsplineOptimizeTraj(ctrl_pts, ts, cost_function, 1, 1);

    t_opt = (ros::Time::now() - t1).toSec();
  }

  if(planning_type_ == 1){
    //benchmark
    t1 = ros::Time::now();
    geo_path_finder_->reset();
    int status = geo_path_finder_->search(start_pt, end_pt, false, 0);
    if (status == Astar::NO_PATH) {
      cout << "[kino replan]: Astar search fail!" << endl;

      // retry searching with discontinuous initial state
      geo_path_finder_->reset();
      status = geo_path_finder_->search(start_pt, end_pt, false, 0);

      if (status == Astar::NO_PATH) {
        cout << "[kino replan]: Can't find path." << endl;
        return false;
      } else {
        cout << "[kino replan]: retry search success." << endl;
      }

    } else {
      cout << "[kino replan]: Astar search success." << endl;
    }

    local_data_.geo_astar_wpts_ = geo_path_finder_->getPath();
    t_search = (ros::Time::now() - t1).toSec();
    t1 = ros::Time::now();
    std::vector<Eigen::Vector3d> astar_grid_path = geo_path_finder_->getPath();
    int temp_index = pp_.max_vel_ * replan_thresh_ * 2.0 / geo_astar_resolution_;
    int astar_size = astar_grid_path.size();
    int forward_n = min(temp_index, astar_size - 1);
    Eigen::Vector3d motion_primitives_goal;


    motion_primitives_goal = astar_grid_path[forward_n];
    // if(forward_n == 0){
    //   motion_primitives_goal = (start_pt + end_pt) / 2;
    // }
    // else{
    //   motion_primitives_goal = astar_grid_path[forward_n];
    // }
    
    motion_primitives_goal = astar_grid_path[forward_n];
    //generate motion primitives
    vector<PolynomialTraj> motion_primitive_buffer;
    motion_primitive_buffer.resize( primitive_num_ * 3); //gnd, aerial_upper, aerial_lower
    double lowest_score = 100000000;
    int lowest_index, lowest_type; // 0:gnd 1:aerial_upper 2:aerial_lower

    double primitive_length = (start_pt.head(2) - motion_primitives_goal.head(2)).norm() / 1.5;
    double vertical_height = 1.2;
    local_data_.Astar_Local_Target_ = motion_primitives_goal;

    vector<Eigen::Vector3d> goal_list;
    Eigen::Vector3d best_goal;
    best_goal.setZero();

    for(int i = 0; i < primitive_num_; i++){

      double mp_orientation = odom_yaw + M_PI / primitive_num_ * i;
      double traj_score = 0;
      double dist = primitive_length;
      double time = 0;
      Eigen::Vector3d goal;
      Eigen::MatrixXd pos(3, 3);
      Eigen::VectorXd t(2);
      //gnd_traj
      PolynomialTraj gnd_traj;
      goal << primitive_length * sin(mp_orientation) + start_pt(0) , primitive_length * cos(mp_orientation) + start_pt(1), start_pt(2);
      goal_list.push_back(goal);
      time = pow(pp_.max_vel_, 2) / pp_.max_acc_ > dist ? sqrt(dist / pp_.max_acc_) : (dist - pow(pp_.max_vel_, 2) / pp_.max_acc_) / pp_.max_vel_ + 2 * pp_.max_vel_ / pp_.max_acc_;
      // time *= 1.2;
      Eigen::Vector3d horizon_dir = ((start_pt - goal).cross(Eigen::Vector3d(0, 0, 1))).normalized();
      Eigen::Vector3d vertical_dir = ((start_pt - goal).cross(horizon_dir)).normalized();
      Eigen::Vector3d random_inserted_pt = (start_pt + goal) / 2 +
                                            (((double)rand()) / RAND_MAX - 0.5) * (start_pt - goal).norm() * goal * 0.8 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989) +
                                            (((double)rand()) / RAND_MAX - 0.5) * (start_pt - goal).norm() * vertical_dir * 0.4 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989);      
      
      random_inserted_pt = (start_pt + goal) / 2 ;

      pos.col(0) = start_pt;
      pos.col(1) = random_inserted_pt;
      pos.col(2) = goal;      
      
      t(0) = t(1) = time / 2;
      
      gnd_traj = PolynomialTraj::minSnapTraj(pos, start_vel, Eigen::Vector3d::Zero(), start_acc, Eigen::Vector3d::Zero(), t);


      if(!geo_path_finder_->checkMotionPrimitive(gnd_traj)) traj_score = 1000000;
      else traj_score = geo_path_finder_->scoreMotionPrimitive(gnd_traj, goal, motion_primitives_goal);
      if(lowest_score > traj_score){
        lowest_score = traj_score;
        lowest_type = 0;
        lowest_index = i;
        best_goal = goal;
      }
      motion_primitive_buffer[i] = gnd_traj;

      //aerial_upper_traj
      PolynomialTraj aerial_upper_traj;
      if(start_pt(2) + vertical_height > motion_primitives_goal[2]){
        goal << primitive_length * sin(mp_orientation) + start_pt(0), primitive_length * cos(mp_orientation) + start_pt(1), motion_primitives_goal[2];
      }        
      else {
        goal << primitive_length * sin(mp_orientation) + start_pt(0), primitive_length * cos(mp_orientation )+ start_pt(1), start_pt(2) + vertical_height;
      }

      goal << primitive_length * sin(mp_orientation) + start_pt(0), primitive_length * cos(mp_orientation )+ start_pt(1), start_pt(2) + vertical_height;
      dist = (start_pt - goal).norm();

      goal_list.push_back(goal);

      time = pow(pp_.max_vel_, 2) / pp_.max_acc_ > dist ? sqrt(dist / pp_.max_acc_) : (dist - pow(pp_.max_vel_, 2) / pp_.max_acc_) / pp_.max_vel_ + 2 * pp_.max_vel_ / pp_.max_acc_;
      // time *= 1.2;
      horizon_dir = ((start_pt - goal).cross(Eigen::Vector3d(0, 0, 1))).normalized();
      vertical_dir = ((start_pt - goal).cross(horizon_dir)).normalized();
      random_inserted_pt = (start_pt + goal) / 2 +
                                            (((double)rand()) / RAND_MAX - 0.5) * (start_pt - goal).norm() * goal * 0.8 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989) +
                                            (((double)rand()) / RAND_MAX - 0.5) * (start_pt - goal).norm() * vertical_dir * 0.4 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989);      

      random_inserted_pt = (start_pt + goal) / 2 ;
      pos.col(0) = start_pt;
      pos.col(1) = random_inserted_pt;
      pos.col(2) = goal;      

      t(0) = t(1) = time / 2;      
      aerial_upper_traj = PolynomialTraj::minSnapTraj(pos, start_vel, Eigen::Vector3d::Zero(), start_acc, Eigen::Vector3d::Zero(), t);
      if(!geo_path_finder_->checkMotionPrimitive(aerial_upper_traj)) traj_score = 1000000;
      else traj_score = geo_path_finder_->scoreMotionPrimitive(aerial_upper_traj, goal, motion_primitives_goal);
      if(lowest_score > traj_score){
        lowest_score = traj_score;
        lowest_type = 1;
        lowest_index = i;
        best_goal = goal;
      }
      motion_primitive_buffer[i + primitive_num_] = aerial_upper_traj;

      //aerial_lower_traj
      PolynomialTraj aerial_lower_traj;
      if(start_pt(2) - vertical_height < 0){
        goal << primitive_length * sin(mp_orientation) + start_pt(0), primitive_length * cos(mp_orientation) + start_pt(1), 0;
      }
      else {
        goal << primitive_length * sin(mp_orientation) + start_pt(0), primitive_length * cos(mp_orientation) + start_pt(1), start_pt(2) - vertical_height;
      }

      dist = (start_pt - goal).norm();
      
      goal_list.push_back(goal);
      time = pow(pp_.max_vel_, 2) / pp_.max_acc_ > dist ? sqrt(dist / pp_.max_acc_) : (dist - pow(pp_.max_vel_, 2) / pp_.max_acc_) / pp_.max_vel_ + 2 * pp_.max_vel_ / pp_.max_acc_;
      // time *= 1.2;
      horizon_dir = ((start_pt - goal).cross(Eigen::Vector3d(0, 0, 1))).normalized();
      vertical_dir = ((start_pt - goal).cross(horizon_dir)).normalized();
      random_inserted_pt = (start_pt + goal) / 2 +
                                            (((double)rand()) / RAND_MAX - 0.5) * (start_pt - goal).norm() * goal * 0.8 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989) +
                                            (((double)rand()) / RAND_MAX - 0.5) * (start_pt - goal).norm() * vertical_dir * 0.4 * (-0.978 / (continous_failures_count_ + 0.989) + 0.989);      
      random_inserted_pt = (start_pt + goal) / 2 ;
      pos.col(0) = start_pt;
      pos.col(1) = random_inserted_pt;
      pos.col(2) = goal;      

      t(0) = t(1) = time / 2;
      aerial_lower_traj = PolynomialTraj::minSnapTraj(pos, start_vel, Eigen::Vector3d::Zero(), start_acc, Eigen::Vector3d::Zero(), t);
      if(!geo_path_finder_->checkMotionPrimitive(aerial_lower_traj)) traj_score = 1000000;  
      else traj_score = geo_path_finder_->scoreMotionPrimitive(aerial_lower_traj, goal, motion_primitives_goal);
      if(lowest_score > traj_score){
        lowest_score = traj_score;
        lowest_type = 2;
        lowest_index = i;
        best_goal = goal;
      }
      motion_primitive_buffer[i + 2 * primitive_num_] = aerial_lower_traj;

    }

    double t_mp = (ros::Time::now() - t1).toSec();
    cout << "t_search: " << t_search << " t_mp:" << t_mp << " best_score: " << lowest_score << endl;
    best_traj = motion_primitive_buffer[lowest_index + lowest_type * primitive_num_];
    //for visualization
    for(int i = 0; i < motion_primitive_buffer.size(); i++){
      vector<Eigen::Vector3d> primitive_wpts;
      primitive_wpts = geo_path_finder_->sampleMotionPrimitive(motion_primitive_buffer[i], 0.001);
      local_data_.motion_primitive_wpts_[i] = primitive_wpts;
    }
    

    local_data_.best_primitive_index_ = lowest_index + lowest_type * primitive_num_;
    local_data_.best_traj_ = best_traj;
    local_data_.duration_ = best_traj.getTimeSum();
    local_data_.compute_time_ = t_search + t_mp;
    continous_failures_count_ = 0;
    
    return true;
  }

  // iterative time adjustment

  t1                    = ros::Time::now();
  NonUniformBspline pos = NonUniformBspline(ctrl_pts, 3, ts);

  double to = pos.getTimeSum();
  pos.setPhysicalLimits(pp_.max_vel_, pp_.max_acc_);
  bool feasible = pos.checkFeasibility(false);

  int iter_num = 0;
  while (!feasible && ros::ok()) {

    feasible = pos.reallocateTime();

    if (++iter_num >= 3) break;
  }

  // pos.checkFeasibility(true);
  // cout << "[Main]: iter num: " << iter_num << endl;

  double tn = pos.getTimeSum();

  cout << "[kino replan]: Reallocate ratio: " << tn / to << endl;
  if (tn / to > 3.0) ROS_ERROR("reallocate error.");

  t_adjust = (ros::Time::now() - t1).toSec();

  // save planned results

  local_data_.position_traj_ = pos;

  double t_total = t_search + t_opt + t_adjust;
  local_data_.compute_time_ = t_total;
  cout << "[kino replan]: time: " << t_total << ", search: " << t_search << ", optimize: " << t_opt
       << ", adjust time:" << t_adjust << endl;

  updateTrajInfo();
  continous_failures_count_ = 0;
  return true;
}

// !SECTION
void FastPlannerManager::updateTrajInfo() {
  local_data_.velocity_traj_     = local_data_.position_traj_.getDerivative();
  local_data_.acceleration_traj_ = local_data_.velocity_traj_.getDerivative();
  local_data_.jerk_traj_         = local_data_.acceleration_traj_.getDerivative();
  local_data_.start_pos_         = local_data_.position_traj_.evaluateDeBoorT(0.0);
  local_data_.duration_          = local_data_.position_traj_.getTimeSum();
  local_data_.traj_id_ += 1;
}

}  // namespace fast_planner
