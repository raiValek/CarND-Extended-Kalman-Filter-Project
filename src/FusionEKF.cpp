#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  noise_ax_ = 9;
  noise_ay_ = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    double px = 0;
    double py = 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      
      px = rho * cos(phi);
      py = rho * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
    }

    MatrixXd H(2,4);
    H <<  1, 0, 0, 0,
          0, 1, 0, 0;

    previous_timestamp_ = measurement_pack.timestamp_;

    VectorXd x(4);
    x << px, py, 0.0, 0.0; 
    MatrixXd P = MatrixXd::Identity(4, 4);
    MatrixXd F = MatrixXd::Identity(4, 4);
    MatrixXd Q = MatrixXd::Identity(4, 4);

    ekf_.Init(x, P, F, H, R_laser_, Q);
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ***************************************************************************/

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;

  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ << 1.0, 0.0, dt, 0.0,
             0.0, 1.0, 0.0, dt,
             0.0, 0.0, 1.0, 0.0,
             0.0, 0.0, 0.0, 1.0;

  double dt_2 = dt*dt;
  double dt_3 = dt_2*dt;
  double dt_4 = dt_3*dt;

  ekf_.Q_ << dt_4/4*noise_ax_, 0.0, dt_3/2*noise_ax_, 0.0,
	     0.0, dt_4/4*noise_ay_, 0, dt_3/2*noise_ay_,
	     dt_3/2*noise_ax_, 0.0, dt_2*noise_ax_, 0.0,
	     0.0, dt_3/2*noise_ay_, 0.0, dt_2*noise_ay_; 

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    double rho = measurement_pack.raw_measurements_[0];
    double phi = measurement_pack.raw_measurements_[1];
    double rho_dot = measurement_pack.raw_measurements_[2];

    VectorXd z(3);
    z << rho, phi, rho_dot;

    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(z);
  }
  else {
    // Laser updates
    double px = measurement_pack.raw_measurements_[0];
    double py = measurement_pack.raw_measurements_[1];

    VectorXd z(2);
    z << px, py;

    ekf_.R_ = R_laser_;

    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
