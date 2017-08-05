#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &P_in, const MatrixXd &F_in,
                        const MatrixXd &H_in, const MatrixXd &R_in, const MatrixXd &Rext_in,
                        const MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Rext_ = Rext_in;
  Q_ = Q_in;

  // creating transposed matrix
  Ht_ = H_.transpose();

  // creating the identity matrix
  long x_size = x_.size();
  I_ = MatrixXd::Identity(x_size, x_size);
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	VectorXd y = z - H_ * x_;

  // support matrices
	MatrixXd S = H_ * P_ * Ht_ + R_;
	MatrixXd K = P_ * Ht_ * S.inverse();

	// new estimate
	x_ = x_ + (K * y);
	P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &Hj) {
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  VectorXd hx = VectorXd(3);

  float sum_sq = sqrt(px * px + py * py);
  hx << sum_sq, atan2(py, px), 0;
  if (sum_sq != 0) {
    hx(2) = (px * vx + py * vy) / sum_sq;
  }

	VectorXd y = z - hx;

  // normalizing polar angle
  while (y(1) < -M_PI) {
    y(1) += 2 * M_PI;
  }
  while (y(1) > M_PI) {
    y(1) -= 2 * M_PI;
  }

  MatrixXd Hjt = Hj.transpose();

	MatrixXd S = Hj * P_ * Hjt + Rext_;
	MatrixXd K = P_ * Hjt * S.inverse();

	// new estimate
	x_ = x_ + (K * y);
	P_ = (I_ - K * Hj) * P_;
}
