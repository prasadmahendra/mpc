#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC
{
public:
  MPC();
  virtual ~MPC();

  void WaypointData2VehicleCoords(const Eigen::VectorXd& current_state, const std::vector<double>& ptsx, const std::vector<double>& ptsy, Eigen::VectorXd& ptsx_conv, Eigen::VectorXd& ptsy_conv);
  void WaypointVectors(const Eigen::VectorXd& ptsx_conv, const Eigen::VectorXd& ptsy_conv, std::vector<double>& next_x_vals, std::vector<double>& next_y_vals);
  void CalculateErrors(const Eigen::VectorXd& current_state);
  vector<double> Solve(const Eigen::VectorXd& current_state);
  
  double ThrottleNext();
  double SteeringAngleNext();
  std::vector<double>& PredictedXVals();
  std::vector<double>& PredictedYVals();

private:
  Eigen::VectorXd polyfit(const Eigen::VectorXd& xvals, const Eigen::VectorXd& yvals, int order);
  double polyeval(const Eigen::VectorXd& coeffs, const double x);
  
  double cte;                               // cross track error; offset from the center of the road
  double epsi;                              // psi error
  double steering_angle_next;               // predicted actuator control; steering angle (δ)
  double throttle_next;                     // predicted actuator control; acceleration (a; throttle + brakes combined in the range of -1 to 1)
  Eigen::VectorXd wayPtPolynomialCoeffs;    // polyfit coeffs of our waypoint coordinates
  std::vector<double> predicted_xvals;      // predicted xvals (vehicles current trajectory)
  std::vector<double> predicted_yvals;      // predicted yvals (vehicles current trajectory)
};

#endif /* MPC_H */
