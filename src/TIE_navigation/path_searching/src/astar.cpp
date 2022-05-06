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



#include <path_searching/astar.h>
#include <sstream>

using namespace std;
using namespace Eigen;

namespace fast_planner {
Astar::~Astar() {
  for (int i = 0; i < allocate_num_; i++) {
    delete path_node_pool_[i];
  }
}

int Astar::search(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt, bool dynamic, double time_start) {
  /* ---------- initialize ---------- */
  cout << "Astar begins!" << endl;
  NodePtr cur_node = path_node_pool_[0];
  cur_node->parent = NULL;
  cur_node->position = start_pt;
  cur_node->index = posToIndex(start_pt);
  cur_node->g_score = 0.0;

  Eigen::Vector3i end_index;
  double time_to_goal;
  end_index = posToIndex(end_pt);
  cur_node->f_score = lambda_heu_ * getEuclHeu(cur_node->position, end_pt);
  cur_node->node_state = IN_OPEN_SET;
  open_set_.push(cur_node);
  use_node_num_ += 1;
  if( cur_node->position[2] < 0 ) cur_node->position[2] = 0;
  if( cur_node->position[2] >= ground_judge_ ){
    cur_node->g_score += aerial_penalty_ * cur_node->position[2] / 2.0;
    cur_node->f_score += aerial_penalty_ * cur_node->position[2] / 2.0;
    cur_node->penalty_g_score = aerial_penalty_ * cur_node->position[2] / 2.0;
    cur_node->motion_state = 1; //flying
  }

  else{
    cur_node->motion_state = 0; //rolling
    cur_node->penalty_g_score = 0;
  } 

  if (dynamic) {
    time_origin_ = time_start;
    cur_node->time = time_start;
    cur_node->time_idx = timeToIndex(time_start);
    expanded_nodes_.insert(cur_node->index, cur_node->time_idx, cur_node);
    // cout << "time start: " << time_start << endl;
  } else
    expanded_nodes_.insert(cur_node->index, cur_node);

  NodePtr neighbor = NULL;
  NodePtr terminate_node = NULL;

  /* ---------- search loop ---------- */
  while (!open_set_.empty()) {
    /* ---------- get lowest f_score node ---------- */
    cur_node = open_set_.top();
    // cout << "pos: " << cur_node->state.head(3).transpose() << endl;
    // cout << "time: " << cur_node->time << endl;
    // cout << "dist: " <<
    // edt_environment_->evaluateCoarseEDT(cur_node->state.head(3),
    // cur_node->time) <<
    // endl;

    /* ---------- determine termination ---------- */

    bool reach_end = abs(cur_node->index(0) - end_index(0)) <= 1 &&
        abs(cur_node->index(1) - end_index(1)) <= 1 && abs(cur_node->index(2) - end_index(2)) <= 1;

    if (reach_end) {
      // cout << "[Astar]:---------------------- " << use_node_num_ << endl;
      // cout << "use node num: " << use_node_num_ << endl;
      // cout << "iter num: " << iter_num_ << endl;
      terminate_node = cur_node;
      retrievePath(terminate_node);
      has_path_ = true;

      return REACH_END;
    }

    /* ---------- pop node and add to close set ---------- */
    open_set_.pop();
    // ROS_WARN_STREAM("g_score: " << cur_node->g_score << "h_score: " << cur_node->f_score - cur_node->g_score);

    cur_node->node_state = IN_CLOSE_SET;
    iter_num_ += 1;

    /* ---------- init neighbor expansion ---------- */

    Eigen::Vector3d cur_pos = cur_node->position;
    Eigen::Vector3d pro_pos;
    double pro_t;

    vector<Eigen::Vector3d> inputs;
    Eigen::Vector3d d_pos;
    // cout << "resolution: " << resolution_ << endl;
    /* ---------- expansion loop ---------- */
    for (double dx = -resolution_; dx <= resolution_ + 1e-3; dx += resolution_)
      for (double dy = -resolution_; dy <= resolution_ + 1e-3; dy += resolution_)
        for (double dz = -resolution_; dz <= resolution_ + 1e-3; dz += resolution_) {
          d_pos << dx, dy, dz;

          if (d_pos.norm() < 1e-3) continue;

          pro_pos = cur_pos + d_pos;

          /* ---------- check if in feasible space ---------- */
          /* inside map range */
          // if (pro_pos(0) <= origin_(0) || pro_pos(0) >= map_size_3d_(0) || pro_pos(1) <= origin_(1) ||
          //     pro_pos(1) >= map_size_3d_(1) || pro_pos(2) <= origin_(2) ||
          //     pro_pos(2) >= map_size_3d_(2)) {
          //   // cout << "outside map" << endl;
          //   continue;
          // }

          /* not in close set */
          Eigen::Vector3i pro_id = posToIndex(pro_pos);
          int pro_t_id = timeToIndex(pro_t);

          NodePtr pro_node =
              dynamic ? expanded_nodes_.find(pro_id, pro_t_id) : expanded_nodes_.find(pro_id);

          if (pro_node != NULL && pro_node->node_state == IN_CLOSE_SET) {
            // cout << "in closeset" << endl;
            continue;
          }

          /* collision free */
          // double dist = dynamic ?
          // edt_environment_->evaluateCoarseEDT(pro_pos, cur_node->time + dt) :
          //                         edt_environment_->evaluateCoarseEDT(pro_pos,
          //                         -1.0);

          if (edt_environment_->sdf_map_->getInflateOccupancy(pro_pos) != 0){
            // cout << "pro_pos: " << pro_pos << endl;
            // cout << "astar grid occ!" << endl;
            continue;
          }

          /* ---------- compute cost ---------- */
          double time_to_goal, tmp_g_score, tmp_f_score, penalty_g_score;
          bool next_motion_state = false;
          tmp_g_score = d_pos.squaredNorm() + cur_node->g_score;
          if(pro_pos[2] > ground_judge_) {
            tmp_g_score -= cur_node->penalty_g_score;
            tmp_g_score += aerial_penalty_ * pro_pos[2] / 2.0;
            penalty_g_score = aerial_penalty_ * pro_pos[2] / 2.0;
            next_motion_state = true;
          }


          tmp_f_score = tmp_g_score + lambda_heu_ * getEuclHeu(pro_pos, end_pt);

          
          if (pro_node == NULL) {
            pro_node = path_node_pool_[use_node_num_];
            pro_node->index = pro_id;
            pro_node->position = pro_pos;
            pro_node->f_score = tmp_f_score;
            pro_node->g_score = tmp_g_score;
            pro_node->parent = cur_node;
            pro_node->node_state = IN_OPEN_SET;
            pro_node->motion_state = next_motion_state;
            pro_node->penalty_g_score = penalty_g_score; 
            if (dynamic) {
              pro_node->time = cur_node->time + 1.0;
              pro_node->time_idx = timeToIndex(pro_node->time);
            }
            open_set_.push(pro_node);

            if (dynamic)
              expanded_nodes_.insert(pro_id, pro_node->time, pro_node);
            else
              expanded_nodes_.insert(pro_id, pro_node);

            use_node_num_ += 1;
            if (use_node_num_ == allocate_num_) {
              cout << "run out of memory." << endl;
              return NO_PATH;
            }
          } else if (pro_node->node_state == IN_OPEN_SET) {
            if (tmp_g_score < pro_node->g_score) {
              // pro_node->index = pro_id;
              pro_node->position = pro_pos;
              pro_node->f_score = tmp_f_score;
              pro_node->g_score = tmp_g_score;
              pro_node->parent = cur_node;
              pro_node->motion_state = next_motion_state;
              pro_node->penalty_g_score = penalty_g_score; 
              if (dynamic) pro_node->time = cur_node->time + 1.0;
            }
          } else {
            cout << "error type in searching: " << pro_node->node_state << endl;
          }

          /* ----------  ---------- */
        }
  }

  /* ---------- open set empty, no path ---------- */
  cout << "open set empty, no path!" << endl;
  cout << "use node num: " << use_node_num_ << endl;
  cout << "iter num: " << iter_num_ << endl;
  return NO_PATH;
}

void Astar::setParam(ros::NodeHandle& nh) {
  nh.param("astar/resolution_astar", resolution_, -1.0);
  nh.param("astar/time_resolution", time_resolution_, -1.0);
  nh.param("astar/lambda_heu", lambda_heu_, -1.0);
  nh.param("astar/margin", margin_, -1.0);
  nh.param("astar/allocate_num", allocate_num_, -1);
  nh.param("astar/weight_goal", weight_goal_, -1.0);
  nh.param("astar/aerial_penalty", aerial_penalty_, 0.0);
  nh.param("astar/ground_judge", ground_judge_, 0.0);
  tie_breaker_ = 1.0 + 1.0 / 10000;

  cout << "margin:" << margin_ << endl;
}

bool Astar::checkMotionPrimitive(PolynomialTraj traj){
  bool result;
  double tau = traj.getTimeSum();
  for (double t = 0.0; t <= tau; t += 0.05){
    Eigen::Vector3d pos, vel, acc;
    pos = traj.evaluate(t);
    if(edt_environment_->evaluateCoarseEDT(pos, -1.0) < 0){
      result = false;
      return result;
    }
  }
  result = true;
  return result;
}
double Astar::scoreMotionPrimitive(PolynomialTraj traj, Eigen::Vector3d start_pos, Eigen::Vector3d goal_pos){
  double nearest_dist = 10000;
  double score;
  bool dist_flag = false;
  double tau = traj.getTimeSum();
  for (double t = 0.0; t <= tau; t += 0.05){
    Eigen::Vector3d pos, vel, acc;
    pos = traj.evaluate(t);
    if(edt_environment_->evaluateCoarseEDT(pos, -1.0) >= margin_ * 1.5 && edt_environment_->evaluateCoarseEDT(pos, -1.0) < 50){
      if(nearest_dist > edt_environment_->evaluateCoarseEDT(pos, -1.0))
      nearest_dist = edt_environment_->evaluateCoarseEDT(pos, -1.0);
      dist_flag = true;
    }
  }
  if(!dist_flag) nearest_dist = 0;
  score = weight_goal_ * (goal_pos - start_pos).norm() - nearest_dist;
  return score;
}
vector<Eigen::Vector3d> Astar::sampleMotionPrimitive(PolynomialTraj traj, double td){
  vector<Eigen::Vector3d> wpts;
  double tau = traj.getTimeSum();
  for (double t = 0.0; t < tau; t += td){
    Eigen::Vector3d pos;
    pos = traj.evaluate(t);
    wpts.push_back(pos);
  }
  return wpts;  
}
void Astar::retrievePath(NodePtr end_node) {
  NodePtr cur_node = end_node;
  path_nodes_.push_back(cur_node);

  while (cur_node->parent != NULL) {
    cur_node = cur_node->parent;
    path_nodes_.push_back(cur_node);
  }

  reverse(path_nodes_.begin(), path_nodes_.end());
}

std::vector<Eigen::Vector3d> Astar::getPath() {
  vector<Eigen::Vector3d> path;
  for (int i = 0; i < path_nodes_.size(); ++i) {
    path.push_back(path_nodes_[i]->position);
  }
  return path;
}

double Astar::getDiagHeu(Eigen::Vector3d x1, Eigen::Vector3d x2) {
  double dx = fabs(x1(0) - x2(0));
  double dy = fabs(x1(1) - x2(1));
  double dz = fabs(x1(2) - x2(2));

  double h;
  int diag = min(min(dx, dy), dz);
  dx -= diag;
  dy -= diag;
  dz -= diag;

  if (dx < 1e-4) {
    h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
  }
  if (dy < 1e-4) {
    h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
  }
  if (dz < 1e-4) {
    h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
  }
  return tie_breaker_ * h;
}

double Astar::getManhHeu(Eigen::Vector3d x1, Eigen::Vector3d x2) {
  double dx = fabs(x1(0) - x2(0));
  double dy = fabs(x1(1) - x2(1));
  double dz = fabs(x1(2) - x2(2));

  return tie_breaker_ * (dx + dy + dz);
}

double Astar::getEuclHeu(Eigen::Vector3d x1, Eigen::Vector3d x2) {
  return tie_breaker_ * (x2 - x1).norm();
}

void Astar::init() {
  /* ---------- map params ---------- */
  this->inv_resolution_ = 1.0 / resolution_;
  inv_time_resolution_ = 1.0 / time_resolution_;
  edt_environment_->getMapRegion(origin_, map_size_3d_);

  cout << "origin_: " << origin_.transpose() << endl;
  cout << "map size: " << map_size_3d_.transpose() << endl;

  /* ---------- pre-allocated node ---------- */
  path_node_pool_.resize(allocate_num_);
  for (int i = 0; i < allocate_num_; i++) {
    path_node_pool_[i] = new Node;
  }

  use_node_num_ = 0;
  iter_num_ = 0;
}

void Astar::setEnvironment(const EDTEnvironment::Ptr& env) {
  this->edt_environment_ = env;
}

void Astar::reset() {
  expanded_nodes_.clear();
  path_nodes_.clear();

  std::priority_queue<NodePtr, std::vector<NodePtr>, NodeComparator0> empty_queue;
  open_set_.swap(empty_queue);

  for (int i = 0; i < use_node_num_; i++) {
    NodePtr node = path_node_pool_[i];
    node->parent = NULL;
    node->node_state = NOT_EXPAND;
  }

  use_node_num_ = 0;
  iter_num_ = 0;
}

std::vector<NodePtr> Astar::getVisitedNodes() {
  vector<NodePtr> visited;
  visited.assign(path_node_pool_.begin(), path_node_pool_.begin() + use_node_num_ - 1);
  return visited;
}

Eigen::Vector3i Astar::posToIndex(Eigen::Vector3d pt) {
  Vector3i idx = ((pt - origin_) * inv_resolution_).array().floor().cast<int>();

  // idx << floor((pt(0) - origin_(0)) * inv_resolution_), floor((pt(1) -
  // origin_(1)) * inv_resolution_),
  //     floor((pt(2) - origin_(2)) * inv_resolution_);

  return idx;
}

int Astar::timeToIndex(double time) {
  int idx = floor((time - time_origin_) * inv_time_resolution_);
}

}  // namespace fast_planner
