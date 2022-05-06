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

#include <plan_manage/kino_replan_fsm.h>

namespace fast_planner {

void KinoReplanFSM::init(ros::NodeHandle& nh) {
  current_wp_  = 0;
  exec_state_  = FSM_EXEC_STATE::INIT;
  have_target_ = false;
  have_odom_   = false;

  /*  fsm param  */
  nh.param("fsm/flight_type", target_type_, -1);
  nh.param("fsm/thresh_replan", replan_thresh_, -1.0);
  nh.param("fsm/thresh_no_replan", no_replan_thresh_, -1.0);

  nh.param("fsm/waypoint_num", waypoint_num_, -1);
  for (int i = 0; i < waypoint_num_; i++) {
    nh.param("fsm/waypoint" + to_string(i) + "_x", waypoints_[i][0], -1.0);
    nh.param("fsm/waypoint" + to_string(i) + "_y", waypoints_[i][1], -1.0);
    nh.param("fsm/waypoint" + to_string(i) + "_z", waypoints_[i][2], -1.0);
  }
  
  /* initialize main modules */
  planner_manager_.reset(new FastPlannerManager);
  planner_manager_->initPlanModules(nh);
  visualization_.reset(new PlanningVisualization(nh));
  odom_yaw_ = 0;
  /* callback */
  exec_timer_   = nh.createTimer(ros::Duration(0.01), &KinoReplanFSM::execFSMCallback, this);
  safety_timer_ = nh.createTimer(ros::Duration(0.01), &KinoReplanFSM::checkCollisionCallback, this);

  waypoint_sub_ =
      nh.subscribe("/waypoint_generator/waypoints", 1, &KinoReplanFSM::waypointCallback, this);
  odom_sub_ = nh.subscribe("/odom_world", 1, &KinoReplanFSM::odometryCallback, this);

  replan_pub_  = nh.advertise<std_msgs::Empty>("/planning/replan", 10);
  new_pub_     = nh.advertise<std_msgs::Empty>("/planning/new", 10);
  bspline_pub_ = nh.advertise<plan_manage::Bspline>("/planning/bspline", 10);
  poly_pub_ = nh.advertise<quadrotor_msgs::PolynomialTrajectory>("/planning/polynomial", 10);

}

void KinoReplanFSM::waypointCallback(const nav_msgs::PathConstPtr& msg) {
  if (msg->poses[0].pose.position.z < -0.1) return;

  cout << "Triggered!" << endl;
  trigger_ = true;

  if (target_type_ == TARGET_TYPE::MANUAL_TARGET) {
    end_pt_ << msg->poses[0].pose.position.x, msg->poses[0].pose.position.y, 0;

  } else if (target_type_ == TARGET_TYPE::PRESET_TARGET) {
    end_pt_(0)  = waypoints_[current_wp_][0];
    end_pt_(1)  = waypoints_[current_wp_][1];
    end_pt_(2)  = waypoints_[current_wp_][2];
    current_wp_ = (current_wp_ + 1) % waypoint_num_;
  }

  visualization_->drawGoal(end_pt_, 0.4, Eigen::Vector4d(0, 1, 0, 1.0));
  end_vel_.setZero();
  have_target_ = true;

  if (exec_state_ == WAIT_TARGET)
    changeFSMExecState(GEN_NEW_TRAJ, "TRIG");
  else if (exec_state_ == EXEC_TRAJ)
    changeFSMExecState(REPLAN_TRAJ, "TRIG");
}

void KinoReplanFSM::odometryCallback(const nav_msgs::OdometryConstPtr& msg) {
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
  odom_yaw_ = ypr(0) * 180 / M_PI;

  have_odom_ = true;
}

void KinoReplanFSM::changeFSMExecState(FSM_EXEC_STATE new_state, string pos_call) {
  string state_str[5] = { "INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "REPLAN_TRAJ", "EXEC_TRAJ" };
  int    pre_s        = int(exec_state_);
  exec_state_         = new_state;
  cout << "[" + pos_call + "]: from " + state_str[pre_s] + " to " + state_str[int(new_state)] << endl;
}

void KinoReplanFSM::printFSMExecState() {
  string state_str[5] = { "INIT", "WAIT_TARGET", "GEN_NEW_TRAJ", "REPLAN_TRAJ", "EXEC_TRAJ" };

  cout << "[FSM]: state: " + state_str[int(exec_state_)] << endl;
}

void KinoReplanFSM::execFSMCallback(const ros::TimerEvent& e) {

  static ros::Time   time_start;

  static int fsm_num = 0;
  fsm_num++;
  if (fsm_num == 100) {
    printFSMExecState();
    if (!have_odom_) cout << "no odom." << endl;
    if (!trigger_) cout << "wait for goal." << endl;
    fsm_num = 0;
  }

  switch (exec_state_) {
    case INIT: {
      if (!have_odom_) {
        return;
      }
      if (!trigger_) {
        return;
      }
      changeFSMExecState(WAIT_TARGET, "FSM");
      break;
    }

    case WAIT_TARGET: {
      if (!have_target_)
        return;
      else {
        changeFSMExecState(GEN_NEW_TRAJ, "FSM");
      }
      break;
    }

    case GEN_NEW_TRAJ: {
      start_pt_  = odom_pos_;
      start_vel_ = odom_vel_;
      start_acc_.setZero();
      if(!time_init_flag){
        time_start = ros::Time::now();
        time_init_flag = true;
      } 
      
      Eigen::Vector3d rot_x = odom_orient_.toRotationMatrix().block(0, 0, 3, 1);

      bool success = callKinodynamicReplan();
      if (success) {
        LocalTrajData* info     = &planner_manager_->local_data_;
        changeFSMExecState(EXEC_TRAJ, "FSM");
      } else {
        // have_target_ = false;
        // changeFSMExecState(WAIT_TARGET, "FSM");
        changeFSMExecState(GEN_NEW_TRAJ, "FSM");
      }
      break;
    }

    case EXEC_TRAJ: {
      /* determine if need to replan */
      LocalTrajData* info     = &planner_manager_->local_data_;
      ros::Time      time_now = ros::Time::now();
      double         t_cur    = (time_now - info->start_time_).toSec();
      double duration;

      Eigen::Vector3d pos, start_pos;
      if(planner_manager_->getPlanningType() != 1){ 
        t_cur                   = min(info->duration_, t_cur);
        pos = info->position_traj_.evaluateDeBoorT(t_cur);
        start_pos = info->start_pos_;
        duration = info->duration_;
      }
      else{
        start_pos = planner_manager_->best_traj.evaluate(0.0);
        duration = planner_manager_->best_traj.getTimeSum();
        if(t_cur >= duration) t_cur = duration;
        pos = planner_manager_->best_traj.evaluate(t_cur - 1e-2);
      }
      /* && (end_pt_ - pos).norm() < 0.5 */

       //finished
       if((end_pt_ - pos).norm() < 0.5){
          //finished
          have_target_ = false;
          changeFSMExecState(WAIT_TARGET, "FSM");
       }
      else if (t_cur > replan_thresh_){
        cout << "fixed time replan!" << endl;
        changeFSMExecState(REPLAN_TRAJ, "FSM");
      }
      else if (t_cur >= duration - 1e-2) {
        if((end_pt_ - pos).norm() < 0.5){
          have_target_ = false;
          changeFSMExecState(WAIT_TARGET, "FSM");
        }
        
        else{
          changeFSMExecState(REPLAN_TRAJ, "FSM");
        }
        return;

      } 
      break;
    }

    case REPLAN_TRAJ: {
      LocalTrajData* info     = &planner_manager_->local_data_;
      ros::Time      time_now = ros::Time::now();
      double         t_cur    = (time_now - info->start_time_).toSec();
      double duration;
      
      Eigen::Vector3d pos;
      if(planner_manager_->getPlanningType() != 1){
        t_cur                   = min(info->duration_, t_cur);
        pos = info->position_traj_.evaluateDeBoorT(t_cur);
        start_pt_  = info->position_traj_.evaluateDeBoorT(t_cur);
        start_vel_ = info->velocity_traj_.evaluateDeBoorT(t_cur);
        start_acc_ = info->acceleration_traj_.evaluateDeBoorT(t_cur);
      }
      else{
        duration = planner_manager_->best_traj.getTimeSum();
        if(t_cur > duration) t_cur = duration;
        start_pt_  = planner_manager_->best_traj.evaluate(t_cur);
        start_vel_ = planner_manager_->best_traj.evaluateVel(t_cur);
        start_acc_ = planner_manager_->best_traj.evaluateAcc(t_cur);
      }


      std_msgs::Empty replan_msg;
      replan_pub_.publish(replan_msg);

      bool success = callKinodynamicReplan();
      if (success) {
        LocalTrajData* info     = &planner_manager_->local_data_;
        changeFSMExecState(EXEC_TRAJ, "FSM");
      } else {
        changeFSMExecState(GEN_NEW_TRAJ, "FSM");
      }
      break;
    }
  }
}

void KinoReplanFSM::checkCollisionCallback(const ros::TimerEvent& e) {
  LocalTrajData* info = &planner_manager_->local_data_;
  double duration;
  if(planner_manager_->getPlanningType() != 1) duration = info->duration_;
  else duration = planner_manager_->best_traj.getTimeSum();;

  if (have_target_) {
    
    auto edt_env = planner_manager_->edt_environment_;

    double dist = planner_manager_->pp_.dynamic_ ?
        edt_env->evaluateCoarseEDT(end_pt_, /* time to program start + */ duration) :
        edt_env->evaluateCoarseEDT(end_pt_, -1.0);

    if (dist <= 0.5) {
      /* try to find a max distance goal around */
      bool            new_goal = false;
      const double    dr = 0.5, dtheta = 30, dz = 0.3;
      double          new_x, new_y, new_z, max_dist = -1.0;
      Eigen::Vector3d goal;

      for (double r = dr; r <= 5 * dr + 1e-3; r += dr) {
        for (double theta = -90; theta <= 270; theta += dtheta) {
          for (double nz = 1 * dz; nz >= -1 * dz; nz -= dz) {

            new_x = end_pt_(0) + r * cos(theta / 57.3);
            new_y = end_pt_(1) + r * sin(theta / 57.3);
            new_z = 0;

            Eigen::Vector3d new_pt(new_x, new_y, new_z);
            dist = planner_manager_->pp_.dynamic_ ?
                edt_env->evaluateCoarseEDT(new_pt, /* time to program start+ */ duration) :
                edt_env->evaluateCoarseEDT(new_pt, -1.0);

            if (dist > max_dist) {
              /* reset end_pt_ */
              goal(0)  = new_x;
              goal(1)  = new_y;
              goal(2)  = new_z;
              max_dist = dist;
            }
          }
        }
      }

      if (max_dist > 0.5) {
        cout << "change goal, replan." << endl;
        end_pt_      = goal;
        have_target_ = true;
        end_vel_.setZero();

        if (exec_state_ == EXEC_TRAJ) {
          changeFSMExecState(REPLAN_TRAJ, "SAFETY");
        }

        visualization_->drawGoal(end_pt_, 0.3, Eigen::Vector4d(1, 0, 0, 1.0));
      } else {
        // have_target_ = false;
        // cout << "Goal near collision, stop." << endl;
        // changeFSMExecState(WAIT_TARGET, "SAFETY");
        cout << "goal near collision, keep retry" << endl;
        changeFSMExecState(REPLAN_TRAJ, "FSM");
        std_msgs::Empty emt;
        replan_pub_.publish(emt);
      }
    }
  }




  if (exec_state_ == FSM_EXEC_STATE::EXEC_TRAJ) {
  /* ---------- check trajectory ---------- */
    double dist;
    int   type = planner_manager_->checkTrajCollision(dist, duration);
    if(type == 2){
      cout << "mission failed!" << endl;
      changeFSMExecState(WAIT_TARGET, "SAFETY");
    }
    if (type == 1) {
      // cout << "current traj in collision." << endl;
      ROS_ERROR("current traj in collision.");
      changeFSMExecState(REPLAN_TRAJ, "SAFETY");
    }
  }
}

bool KinoReplanFSM::callKinodynamicReplan() {
  bool plan_success =
      planner_manager_->kinodynamicReplan(start_pt_, start_vel_, start_acc_, end_pt_, end_vel_);
  
  auto info = &planner_manager_->local_data_;
  if (plan_success && (planner_manager_->getPlanningType() != 1)) {

    /* publish traj */
    plan_manage::Bspline bspline;
    bspline.order      = 3;
    bspline.start_time = info->start_time_;
    bspline.traj_id    = info->traj_id_;

    Eigen::MatrixXd pos_pts = info->position_traj_.getControlPoint();

    for (int i = 0; i < pos_pts.rows(); ++i) {
      geometry_msgs::Point pt;
      pt.x = pos_pts(i, 0);
      pt.y = pos_pts(i, 1);
      pt.z = pos_pts(i, 2);
      bspline.pos_pts.push_back(pt);
    }

    Eigen::VectorXd knots = info->position_traj_.getKnot();
    for (int i = 0; i < knots.rows(); ++i) {
      bspline.knots.push_back(knots(i));
    }

    bspline_pub_.publish(bspline);

    /* visulization */
    auto plan_data = &planner_manager_->plan_data_;
    visualization_->drawGeometricPath(plan_data->kino_path_, 0.075, Eigen::Vector4d(1, 1, 0, 0.4));
    visualization_->drawBspline(info->position_traj_, 0.1, Eigen::Vector4d(1.0, 0, 0.0, 1), true, 0.2,
                                Eigen::Vector4d(0, 0, 1, 1));                               
    return true;
  } 
  else if(plan_success && (planner_manager_->getPlanningType() == 1)){
    polynomialTrajConverter(planner_manager_->best_traj);
    if(info-> geo_astar_wpts_.size() > 0) visualization_->drawGeometricPath(info->geo_astar_wpts_, 0.2, Eigen::Vector4d(0, 0, 1, 1.0)); 
    visualization_->drawGoal(info->Astar_Local_Target_, 0.5, Eigen::Vector4d(0, 0, 0, 1.0));
    visualization_->drawMotionPrimitive(info->motion_primitive_wpts_, info->best_primitive_index_, 
                                0.05, Eigen::Vector4d(0, 1, 0, 1.0), Eigen::Vector4d(1, 0, 0, 1.0));
    visualization_->drawBspline(info->position_traj_, 0.1, Eigen::Vector4d(1.0, 0, 0.0, 1), true, 0.2, Eigen::Vector4d(1, 0, 0, 1));
    // cout << "visualization finished" << endl;
    return true;
  }
  else {
    cout << "generate new traj fail." << endl;
    return false;
  }
}
void KinoReplanFSM::polynomialTrajConverter(PolynomialTraj &traj)
{
    quadrotor_msgs::PolynomialTrajectory trajMsg;
    trajMsg.header.stamp = ros::Time::now();
    static uint32_t traj_id = 0;
    traj_id++;
    trajMsg.trajectory_id = traj_id;
    trajMsg.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;
    trajMsg.num_order = 5;
    trajMsg.num_segment = 2;
    trajMsg.header.frame_id = "world";
    Eigen::Vector3d initialVel, finalVel;
    trajMsg.start_yaw = 0.0;
    trajMsg.final_yaw = 0.0;
    vector<double> times = traj.getTimes();
    trajMsg.time = times;
    vector<vector<double>> x_coef = traj.getCoef(0);
    vector<vector<double>> y_coef = traj.getCoef(1);
    vector<vector<double>> z_coef = traj.getCoef(2);
    for (int i = 0; i < 2; ++i)
    {

      trajMsg.order.push_back(5);
      for (size_t j = 0; j <= trajMsg.order[i]; ++j)
      {
        trajMsg.coef_x.push_back(x_coef[i][j]);
        trajMsg.coef_y.push_back(y_coef[i][j]);
        trajMsg.coef_z.push_back(z_coef[i][j]);
      }
    }
    
    trajMsg.mag_coeff = 1.0;
    trajMsg.debug_info = "";
    poly_pub_.publish(trajMsg);
}
// KinoReplanFSM::
}  // namespace fast_planner
