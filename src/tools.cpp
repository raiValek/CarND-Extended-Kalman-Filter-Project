#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd rmse(4);
  rmse << 0.0,0.0,0.0,0.0;

  if(estimations.size() != ground_truth.size()) {
    cout << "CalculateRMSE: Size of estimation data and ground truth data differ" << endl;
    return rmse;
  }

  if(estimations.empty()) {
    cout << "CalculateRMSE: Estimation data is empty" << endl;
    return rmse;
  }


  for(size_t i=0; i < estimations.size(); ++i) {
    VectorXd dif = estimations[i]-ground_truth[i];
    rmse += (VectorXd)(dif.array()*dif.array());
  }

  rmse = rmse.array() / estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;
}

VectorXd Tools::CalculateFloatingRMSE(const VectorXd &estimation,
                                      const VectorXd &ground_truth) {
  VectorXd dif = estimation-ground_truth;
  sum_squares_ += (VectorXd)(dif.array()*dif.array());
  ++num_estimations_;
  return (sum_squares_.array() / num_estimations_).sqrt();
}

MatrixXd Tools::CalculateJacobian(const Vector4d& x_state) {

  MatrixXd Hj(3,4);
  Hj << 0,0,0,0,
        0,0,0,0,
        0,0,0,0;

  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  if(px == 0.0 && px == 0.0){
    cout << "CalculateJacobian - Error - Division by Zero" << endl;
    return Hj;
  }
  
  double px_py = pow(px,2)+pow(py,2);
  double sqr_px_py = sqrt(px_py);
  double sqr_px_py_3 = sqrt(pow(px_py,3));

  Hj << px/sqr_px_py, py/sqr_px_py, 0.0, 0.0,
       -py/px_py, px/px_py, 0, 0,
        py*(vx*py-vy*px)/sqr_px_py_3,
        px*(vy*px-vx*py)/sqr_px_py_3,
        px/sqr_px_py, py/sqr_px_py;
}
