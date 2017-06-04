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
  vector<double> Solve(const Eigen::VectorXd& state);
  
  double ThrottleNext();
  double SteeringAngleNext();
  std::vector<double>& PredictedXVals();
  std::vector<double>& PredictedYVals();

  double cte;       // cross track error; offset from the center of the road
  double epsi;      // psi error

private:
  Eigen::VectorXd polyfit(const Eigen::VectorXd& xvals, const Eigen::VectorXd& yvals, int order);
  double polyeval(const Eigen::VectorXd& coeffs, const double x);
  
  Eigen::VectorXd wayPtPolynomialCoeffs;
  double steering_angle_next;
  double throttle_next;
  std::vector<double> predicted_xvals;
  std::vector<double> predicted_yvals;
};

#endif /* MPC_H */
