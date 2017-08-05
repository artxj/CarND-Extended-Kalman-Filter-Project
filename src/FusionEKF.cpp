#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

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
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  tools = Tools();
  ekf_ = KalmanFilter();

  // state covariance matrix
  MatrixXd P = MatrixXd(4, 4);
	P << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

  // laser measurement matrix
  MatrixXd H = MatrixXd(2, 4);
	H << 1, 0, 0, 0,
			  0, 1, 0, 0;

	// the initial transition matrix
	MatrixXd F = MatrixXd(4, 4);
	F << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

  // the initial process covariance matrix
  MatrixXd Q = MatrixXd(4, 4);
  Q << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

  ekf_.Init(VectorXd(4), P, F, H, R_laser_, R_radar_, Q);

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
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */

    switch (measurement_pack.sensor_type_) {
      case MeasurementPackage::RADAR: {
        float rho = measurement_pack.raw_measurements_[0];
        float phi = measurement_pack.raw_measurements_[1];
        float px, py;
        tools.CalculateCartesian(rho, phi, px, py);
        ekf_.x_ << px, py, 0, 0;
        break;
      }
      case MeasurementPackage::LASER: {
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        break;
      }
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
 	previous_timestamp_ = measurement_pack.timestamp_;
  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // updating state transition matrix F
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // updating the covariance matrix
  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;
  float dt3_div2 = dt3 / 2;
  float dt4_div4 = dt4 / 4;

  ekf_.Q_(0, 0) = noise_ax_ * dt4_div4;
  ekf_.Q_(0, 2) = noise_ax_ * dt3_div2;
  ekf_.Q_(1, 1) = noise_ay_ * dt4_div4;
  ekf_.Q_(1, 3) = noise_ay_ * dt3_div2;
  ekf_.Q_(2, 0) = noise_ax_ * dt3_div2;
  ekf_.Q_(2, 2) = noise_ax_ * dt2;
  ekf_.Q_(3, 1) = noise_ay_ * dt3_div2;
  ekf_.Q_(3, 3) = noise_ay_ * dt2;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  switch (measurement_pack.sensor_type_) {
    case MeasurementPackage::RADAR: {
      VectorXd z = VectorXd(3);
      z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],
                      measurement_pack.raw_measurements_[2];
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.UpdateEKF(z, Hj_);
      break;
    }
    case MeasurementPackage::LASER: {
      VectorXd z = VectorXd(2);
      z << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1];
      ekf_.Update(z);
      break;
    }
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
