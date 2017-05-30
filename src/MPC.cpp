#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"

using CppAD::AD;

/*
 MPC algorithm:
 Setup:
 
 Define the length of the trajectory, N, and duration of each timestep, dt.
 Define vehicle dynamics and actuator limitations along with other constraints.
 Define the cost function.
 Loop:
 
 We pass the current state as the initial state to the model predictive controller.
 We call the optimization solver. Given the initial state, the solver will return the vector of control inputs that minimizes the cost function. The solver we'll use is called Ipopt.
 We apply the first control input to the vehicle.
 Back to 1.
 
 Vehicle state:
 
 x          = vehicle x position
 y          = vehicle y position
 ψ (psi)    = yaw angle or orientation
 v          = velocity
 δ          = actuator controls; steering angle
 a          = actuator controls; acceleration (throttle + brakes combined in the range of -1 to 1)
 
 State changes over time (dt)
 
 x and y:
 ========
 
 x(t+dt) = x(t) + v(t) * cos(ψ) * dt
 y(t+dt) = y(t) + v(t) * sin(ψ) * dt
 
 orientation/yaw (psi):
 =====================
 
 ψ(t+dt) = ψ(t) + ​v(t)/L​f * δ * dt
 
 In a nutshell, we add a multiplicative factor of the steering angle, δ to ψ. L​f measures the distance between the front of the vehicle and its center of gravity. 
 The larger the vehicle, the slower the turn rate (and hence the inverse relationship). v (velocity) is a factor here since at higher speeds you turn quicker than at lower speeds
 
 velocity (v):
 ============
 v(t+dt) = v(t) + a(t) * dt
 
 */


// Set the timestep length and duration
size_t N = 10;
double dt = 0.05;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// The solver takes all the state variables and actuator variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
int px_range_begin       = 0;
int py_range_begin       = px_range_begin    + N;
int psi_range_begin      = py_range_begin    + N;
int v_range_begin        = psi_range_begin   + N;
int cte_range_begin      = v_range_begin     + N;
int epsi_range_begin     = cte_range_begin   + N;
int steering_range_begin = epsi_range_begin  + N;
int throttle_range_begin = steering_range_begin + (N - 1);


// Multipliers for the cost computation
const double cost_state_cte    = 300.0;
const double cost_state_epsi   = 50.0;
const double cost_state_v      = 1.0;
const double cost_val_steering = 200.0;
const double cost_val_throttle = 50.0;
const double cost_seq_steering = 5000.0;
const double cost_seq_throttle = 100.0;

class FG_eval {
 public:
  Eigen::VectorXd coeffs; // Fitted polynomial coefficients
  FG_eval(Eigen::VectorXd coeffs)
  {
    this->coeffs = coeffs;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  
  void operator()(ADvector& fg, const ADvector& vars)
  {
    // MPC Implementation ...

    std::cout << "fg: " << fg.size() << std::endl;
    std::cout << "vars: " << vars.size() << std::endl;
    
    // `fg` is a vector containing the cost and constraints.
    // `vars` is a vector containing the variable values (state & actuators).
    
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    
    // fg[0] is the cost (error) that we would like to minimize. For example, measuring the offset from the center of the lane,
    // where the center of the lane can be called the reference, or desired, state.
    // cte and eψ (epsi):
    //    Ideally, both of these errors would be 0 - there would be no difference from the actual vehicle position and heading to the desired position and heading.
    //    So our desired cte and desired epsi values are 0.0
    
    double desired_cte = 0.0;
    double desired_epsi = 0.0;

    // Our goal is to move the vehicle from A to B, then coming to a halt in the middle of the reference trajectory is a big problem!
    // A simple solution is to capture the velocity error in the cost function. This will penalize the vehicle for not maintaining the reference velocity.
    // Hence we will set our desired velocity to a speed that we expect our car to maintain throughout the track (Another option is to measure the euclidean distance
    // between the current position of the vehicle and the destination and adding that to the cost.)
    
    double desired_velocity = 40;   // go speed racer, go!
    
    fg[0] = 0.0;
    for (int t = 0; t < N; t++)
    {
      fg[0] += cost_state_cte * pow(vars[cte_range_begin + t]  - desired_cte,  2);
      fg[0] += cost_state_epsi * pow(vars[epsi_range_begin + t] - desired_epsi, 2);
      fg[0] += cost_state_v   * pow(vars[v_range_begin + t]    - desired_velocity, 2);
    }
    
    // The cost function is not limited to the state, we could also include the control input! The reason we would do this is to allow us to penalize the magnitude of
    // the input as well as the change-rate. If we want to change lanes, for example, we would have a large cross-track error, but we wouldn't want to jerk the steering
    // wheel as hard as we can. We could add the control input magnitude like this:

    for (int t = 0; t < N - 1; t++)
    {
      fg[0] += cost_val_steering * pow(vars[steering_range_begin + t], 2);
      fg[0] += cost_val_throttle * pow(vars[throttle_range_begin + t], 2);
    }
    
    // We still need to capture the change-rate of the control input to add some temporal smoothness.
    // This additional term in the cost function captures the difference between the next actuator state and the current one:
    
    for (int t = 0; t < N - 2; t++)
    {
      fg[0] += cost_seq_steering * pow(vars[steering_range_begin + t], 2);
      fg[0] += cost_seq_throttle * pow(vars[throttle_range_begin + t], 2);
    }

    // Initial constraints
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`. This bumps up the position of all the other values.
    fg[1 + px_range_begin] = vars[px_range_begin];
    fg[1 + py_range_begin] = vars[py_range_begin];
    fg[1 + psi_range_begin] = vars[psi_range_begin];
    fg[1 + v_range_begin] = vars[v_range_begin];
    fg[1 + cte_range_begin] = vars[cte_range_begin];
    fg[1 + epsi_range_begin] = vars[epsi_range_begin];
    
    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) {
      // The state at time t+1 .
      AD<double> x1 = vars[px_range_begin + i + 1];
      AD<double> y1 = vars[py_range_begin + i + 1];
      AD<double> psi1 = vars[psi_range_begin + i + 1];
      AD<double> v1 = vars[v_range_begin + i + 1];
      AD<double> cte1 = vars[cte_range_begin + i + 1];
      AD<double> epsi1 = vars[epsi_range_begin + i + 1];
      
      // The state at time t.
      AD<double> x0 = vars[px_range_begin + i];
      AD<double> y0 = vars[py_range_begin + i];
      AD<double> psi0 = vars[psi_range_begin + i];
      AD<double> v0 = vars[v_range_begin + i];
      AD<double> cte0 = vars[cte_range_begin + i];
      AD<double> epsi0 = vars[epsi_range_begin + i];
      
      // Only consider the actuation at time t.
      AD<double> delta0 = vars[steering_range_begin + i];
      AD<double> a0 = vars[throttle_range_begin + i];
      
      AD<double> f0 = coeffs[0] + coeffs[1] * x0;
      AD<double> psides0 = CppAD::atan(coeffs[1]);
      
      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[2 + px_range_begin + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[2 + py_range_begin + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[2 + psi_range_begin + i] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[2 + v_range_begin + i] = v1 - (v0 + a0 * dt);
      fg[2 + cte_range_begin + i] =
      cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[2 + epsi_range_begin + i] =
      epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

// Solve the model given an initial state and polynomial coefficients.
// Return the first actuatotions.
vector<double> MPC::Solve(Eigen::VectorXd current_state) {
  double px     = current_state[0];
  double py     = current_state[1];
  double ψ      = current_state[2];
  double v      = current_state[3];
  double δ      = current_state[4];
  double a      = current_state[5];
  
  int state_size = 6;
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // n_vars: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  // 4 * 10 + 2 * 9

  size_t n_vars = state_size * N + 2 * (N-1);

  // Set the number of constraints
  
  size_t n_constraints = state_size * N;

  // Initial value of the independent variables. SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++)
  {
    vars[i] = 0;
  }
  
  // initial state ...
  // vars is a vehicle state (state + actuators) vector representation for N timesteps
  
  vars[px_range_begin]        = px;
  vars[py_range_begin]        = py;
  vars[psi_range_begin]       = ψ;
  vars[v_range_begin]         = v;
  vars[cte_range_begin]       = cte;
  vars[epsi_range_begin]      = epsi;
  vars[steering_range_begin]     = 0.0;
  vars[throttle_range_begin]  = 0.0;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set lower and upper limits for variables.
  // State values *except* actuator controls don't have an expected limit. Set these to double::min and max values
  // Set px, py, ψ, v, cte, epsi upper and lower bound values to double::min and double::max values respectively.
  for (int i = 0; i < steering_range_begin; i++)
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  
  // steering upper and lower bounds ...
  // The upper and lower limits of steering angles are set to -25 and 25
  // degrees (values in radians).
  
  for (int i = steering_range_begin; i < throttle_range_begin; i++)
  {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // throttle upper and lower bounds ...
  for (int i = throttle_range_begin; i < n_vars; i++)
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = +1.0;
  }
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++)
  {
    constraints_lowerbound[i] = 0.0;
    constraints_upperbound[i] = 0.0;
  }
  
  // ... initial state constraints.
  constraints_lowerbound[px_range_begin]    = px;
  constraints_upperbound[px_range_begin]    = px;

  constraints_lowerbound[py_range_begin]    = py;
  constraints_upperbound[py_range_begin]    = py;

  constraints_lowerbound[psi_range_begin]  = ψ;
  constraints_upperbound[psi_range_begin]  = ψ;

  constraints_lowerbound[v_range_begin]    = v;
  constraints_upperbound[v_range_begin]    = v;

  constraints_lowerbound[cte_range_begin]  = cte;
  constraints_upperbound[cte_range_begin]  = cte;

  constraints_lowerbound[epsi_range_begin] = epsi;
  constraints_upperbound[epsi_range_begin] = epsi;

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // object that computes objective and constraints
  FG_eval fg_eval(wayPtPolynomialCoeffs);
  
  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound, constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Solved (predicted) waypoint ...
  predicted_xvals = vector<double>();
  predicted_yvals = vector<double>();
  
  for (int i = 0; i < N; i++)
  {
    predicted_xvals.push_back(solution.x[px_range_begin + i]);
    predicted_yvals.push_back(solution.x[py_range_begin + i]);
  }
  
  // Return the first actuator values. The variables can be accessed with `solution.x[i]`.
  steering_angle_next = -1.0 * solution.x[steering_range_begin];
  throttle_next       = solution.x[throttle_range_begin];
  
  return { solution.x[px_range_begin + 1], solution.x[py_range_begin + 1],
           solution.x[psi_range_begin + 1], solution.x[v_range_begin + 1],
           solution.x[cte_range_begin + 1], solution.x[epsi_range_begin + 1],
           solution.x[steering_range_begin], solution.x[throttle_range_begin]};
}

void MPC::WaypointData2VehicleCoords(const Eigen::VectorXd& current_state, const std::vector<double>& ptsx, const std::vector<double>& ptsy, Eigen::VectorXd& ptsx_conv, Eigen::VectorXd& ptsy_conv)
{
  assert(ptsx.size() == ptsy.size() && ptsx_conv.size() == ptsy_conv.size());
  
  auto px  = current_state[0];
  auto py  = current_state[1];
  auto psi = current_state[2];
  
  int pts_max = ptsx.size();
  for(int i = 0; i < pts_max; i++) {
    auto dx = ptsx[i] - px;
    auto dy = ptsy[i] - py;
    
    ptsx_conv[i] = dx * cos(-psi) - dy * sin(-psi);
    ptsy_conv[i] = dy * cos(-psi) + dx * sin(-psi);
  }
}

void MPC::CalculateErrors(Eigen::VectorXd current_state) {
  double px = current_state[0];
  double py = current_state[1];
  
  // errors ...
  // cross track error or CTE:
  // CTE is y (read: y is the offset from the center of the road) value calculated at car coord X value == 0.0
  // Our polynomial func f = wayPtPolynomialCoeffs[3] * pow(x, 3) + wayPtPolynomialCoeffs[2] * pow(x, 2) + wayPtPolynomialCoeffs[1] * pow(x, 1) + wayPtPolynomialCoeffs[0];
  
  cte = polyeval(wayPtPolynomialCoeffs, 0.0) - py;
  
  // orientation error
  // eψ(​t) = ψ(​t) − ψdes(​t)
  // where:
  //   eψ(​t)    = orientation error at time (t)
  //   ψ(​t)     = orientation at time (t)
  //   ψdes(​t)  = desired orientation at time (t)
  // We already know ψ(​t)​, because it’s part of our state. We don’t yet know ψdes(​t) (desired psi) - all we have so far is a polynomial to follow. ψdes(​t)
  // can be calculated as the tangential angle of the polynomial f evaluated at x(​t), arctan(f′(x​t)). f​′ is the derivative of the polynomial.
  // f' is = 3 * wayPtPolynomialCoeffs[3] * pow(x, 2) + 2 * wayPtPolynomialCoeffs[2] * pow(x, 1) * wayPtPolynomialCoeffs[1]
  // therefore ...
  
  epsi = -atan(wayPtPolynomialCoeffs[1]);
}

void MPC::WaypointVectors(const Eigen::VectorXd& ptsx_conv, const Eigen::VectorXd& ptsy_conv, std::vector<double>& next_x_vals, std::vector<double>& next_y_vals) {
  assert(ptsx_conv.size() == ptsy_conv.size() && ptsx_conv.size() == ptsy_conv.size());
  
  // fit a third order polynomial to waypoints (converted) vectors & get the polynomial coefficients ...
  wayPtPolynomialCoeffs = polyfit(ptsx_conv, ptsy_conv, 3);
  for (int i = 0; i < next_x_vals.size(); i++) {
    next_x_vals[i] = 5.0 * (double)i;
    next_y_vals[i] = polyeval(wayPtPolynomialCoeffs, next_x_vals[i]);
  }
}

// Fit a polynomial.
// Adapted from https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd MPC::polyfit(const Eigen::VectorXd& xvals, const Eigen::VectorXd& yvals, int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);
  
  for (int i = 0; i < xvals.size(); ++i) {
    A(i, 0) = 1.0;
  }
  
  for (int j = 0; j < xvals.size(); ++j) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }
  
  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

double MPC::polyeval(const Eigen::VectorXd& coeffs, const double x)
{
  double y = 0.0;
  
  for (int i=0; i < coeffs.size(); i++)
  {
    y += coeffs[i] * pow(x, i);
  }
  
  return y;
}

double MPC::ThrottleNext() {
  return throttle_next;
}

double MPC::SteeringAngleNext() {
  return steering_angle_next;
}

std::vector<double>& MPC::PredictedXVals() {
  return predicted_xvals;
}

std::vector<double>& MPC::PredictedYVals() {
  return predicted_yvals;
}


