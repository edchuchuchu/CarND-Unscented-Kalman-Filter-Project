#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if (estimations.size() == 0){
      cout << "Error: no estimations vector!" << endl;
      return rmse;
  }
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()){
      cout << "Error: Vector size mismatch!" << endl;
      return rmse;
  }   

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
      // ... your code here
      VectorXd diff = (estimations[i] - ground_truth[i]).array().pow(2);
      rmse += diff;
  }

  //calculate the mean
  rmse /= estimations.size();
  //calculate the squared root
  rmse = rmse.array().sqrt();
  //return the result
  return rmse;   
}
