#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

class KalmanFilter {
public:

  // state vector
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transition matrix
  Eigen::MatrixXd F_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_;

  // measurement covariance matrix for EKF
  Eigen::MatrixXd Rext_;

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_in Measurement covariance matrix for standard KF
   * @param Rext_in Measurement covariance matrix for extended KF
   * @param Rext_in Process covariance matrix
   */
  void Init(const Eigen::VectorXd &x_in, const Eigen::MatrixXd &P_in, const Eigen::MatrixXd &F_in,
      const Eigen::MatrixXd &H_in, const Eigen::MatrixXd &R_in, const Eigen::MatrixXd &Rext_in,
      const Eigen::MatrixXd &Q_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   * @param Hj Jacobian matrix for current step
   */
  void UpdateEKF(const Eigen::VectorXd &z, const Eigen::MatrixXd &Hj);

private:
  // Identity matrix
  Eigen::MatrixXd I_;

  // Transposed H matrix
  Eigen::MatrixXd Ht_;
  
};

#endif /* KALMAN_FILTER_H_ */
