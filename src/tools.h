#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector4d;
using namespace std;

/**
* helper class for convenience
*/
class Tools {
public:
  Tools(){
    InitRMSE();
  }

  /**
  * Initialize RMSE Calculator. Call before starting a new data set.
  */
  inline void InitRMSE(){
    sum_squares_ = VectorXd(4);
    sum_squares_ << 0.0, 0.0, 0.0, 0.0;
    num_estimations_ = 0;
  }

  /**
  * A helper method to calculate RMSE.
  */
  static VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate RMSE for one datapoint at a time.
  */
  VectorXd CalculateFloatingRMSE(const VectorXd &estimation, const VectorXd &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static MatrixXd CalculateJacobian(const Vector4d& x_state);

private:
  VectorXd sum_squares_;
  long num_estimations_;
};

#endif /* TOOLS_H_ */
