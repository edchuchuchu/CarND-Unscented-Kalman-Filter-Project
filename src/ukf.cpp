#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	is_initialized_ = false;

	previous_timestamp_ = 0;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; //30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;  
  
  //initial matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //initial vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);  
  
  // initial NIS for radar
  NIS_radar_ = 0;

  // initial NIS for laser
  NIS_laser_ = 0;  
  
	//state covariance matrix P  
  P_ <<    1,   0,    0,   0,   0,
        0,   1,    0,   0,   0,
        0,   0,    1,   0,   0,
        0,   0,    0,   1,   0,
        0,   0,    0,   0,   1;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  
  if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float px = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
      float py = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
      x_ << px, py, 0 , 0, 0;
    }
    else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }  
 	//compute the time elapsed between the current and previous measurements 
	double delta_t_ = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;  
  Prediction(delta_t_);
  // print the output
  // cout << "x_ = " << x_ << endl;
  // cout << "P_ = " << P_ << endl;  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Laser updates
    UpdateLidar(meas_package);
  } 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //Step1. Generate Sigma Points 
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.fill(0);
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) << pow(std_a_, 2), 0 , 0, pow(std_yawdd_, 2);
  
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++ ){
      Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }  
  
  //Step2. Predict Sigma Points
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
      double px = Xsig_aug(0, i);
      double py = Xsig_aug(1, i);
      double v = Xsig_aug(2, i);
      double yaw = Xsig_aug(3, i);
      double yawd = Xsig_aug(4, i);
      double v_a = Xsig_aug(5, i);
      double v_yawd = Xsig_aug(6, i);
      
      //avoid division by zero
      if(yawd != 0){
          Xsig_pred_(0, i) = (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw));
          Xsig_pred_(1, i) = (v / yawd) * (-cos(yaw + yawd * delta_t) + cos(yaw));
      }else{
          Xsig_pred_(0, i) = v  * cos(yaw) * delta_t;
          Xsig_pred_(1, i) = v  * sin(yaw) * delta_t;          
      }
      //write predicted sigma points into right column
      Xsig_pred_(0, i) += px + 0.5 * pow(delta_t, 2) * cos(yaw) * v_a;
      Xsig_pred_(1, i) += py + 0.5 * pow(delta_t, 2) * sin(yaw) * v_a;
      Xsig_pred_(2, i) = v + delta_t * v_a;
      Xsig_pred_(3, i) = yaw + delta_t * yawd + 0.5 * pow(delta_t, 2) * v_yawd;
      Xsig_pred_(4, i) = yawd + delta_t * v_yawd ;

  } 
  //Step3. Oreduct Mean and Covariance
  //set weights
  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
      if (i==0){
        weights_(i) = lambda_/ (lambda_ + n_aug_);  
      }else{
        weights_(i) = 1 / (2 * (lambda_ + n_aug_));  
      }
      x_ += weights_(i) * Xsig_pred_.col(i);
  }
  P_.fill(0.0);
  //predict state covariance matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++){

  // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;  
    
    P_ += weights_(i) * x_diff * x_diff.transpose() ;
  }  
  // print the output
  // cout << "x_ = " << x_ << endl;
  // cout << "P_ = " << P_ << endl; 
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2; 
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  //transform sigma points into measurement space
  //calculate mean predicted measurement
  z_pred.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);

    Zsig.col(i) << px, py;  

    z_pred += weights_(i) * Zsig.col(i);
  }
  //calculate measurement covariance matrix S   
  S.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI; 
    
    S += weights_(i) * z_diff * z_diff.transpose();
  }  
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0);
  R(0,0) = pow(std_laspx_, 2);
  R(1,1) = pow(std_laspy_, 2);
  S += R;  
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;  
  MatrixXd K = Tc * S.inverse();
  VectorXd z = meas_package.raw_measurements_;
  
  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();   
  // Calculate Normalized Innovation Squared (NIS)
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3; 
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  //transform sigma points into measurement space
  //calculate mean predicted measurement
  z_pred.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    if ((px!= 0)&&(py!= 0)){
      Zsig.col(i) << sqrt(pow(px, 2) + pow(py, 2)), atan2(py, px), (px * cos(yaw) * v + py * sin(yaw) * v)/(sqrt(pow(px, 2) + pow(py, 2)));  
    }else{
      Zsig.col(i) << 0, 0, 0;  
    }
    z_pred += weights_(i) * Zsig.col(i);
  }
  //calculate measurement covariance matrix S   
  S.fill(0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI; 
    
    S += weights_(i) * z_diff * z_diff.transpose();
  }  
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0);
  R(0,0) = pow(std_radr_, 2);
  R(1,1) = pow(std_radphi_, 2);
  R(2,2) = pow(std_radrd_, 2);
  S += R;  
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  //calculate cross correlation matrix
  Tc.fill(0);
  for (int i=0; i < 2 * n_aug_ + 1; i++) {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;  
  MatrixXd K = Tc * S.inverse();
  VectorXd z = meas_package.raw_measurements_;
  
  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();   
  // Calculate Normalized Innovation Squared (NIS)
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff; 
}
