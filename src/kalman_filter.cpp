#include "kalman_filter.h" 

#define PI 3.14159265358979323846

using namespace Eigen;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Ht_ = H_.transpose();
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * Ht_ + R_;
  MatrixXd K = P_ * Ht_ * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  size_t x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  MatrixXd Hj = Tools::CalculateJacobian(x_);

  // Transformation of cartesian state vector to polar coordinates
  Vector3d hx;

  hx(0) = sqrt(pow(px,2) + pow(py,2));
  
  if(px == 0.0){
    hx(1) = PI/2;
  }
  else{
    hx(1) = atan2(py, px); 
  }
  
  if(hx(0) < 0.000001){
    hx(2) = 0.0;
  }
  else{
    hx(2) = (px*vx+py*vy)/hx(0);
  }

  // Update
  MatrixXd Hjt = Hj.transpose();

  VectorXd y = z - hx;

  while(y(1) > PI) y(1) -= 2*PI;
  while(y(1) < -PI) y(1) += 2*PI;

  MatrixXd S = Hj * P_ * Hjt + R_;
  MatrixXd K = P_ * Hjt * S.inverse();
  
//new estimate
  x_ = x_ + (K * y);
  size_t x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
