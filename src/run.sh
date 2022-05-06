roslaunch plan_manage kino_replan.launch & sleep 1;
roslaunch poly_traj_server traj_server.launch & sleep 1;
wait;