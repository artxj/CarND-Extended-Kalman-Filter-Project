#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

  int size = estimations.size();
  assert(size != 0 && ground_truth.size() == size);

	for (int i = 0; i < size; ++i) {
    VectorXd res = estimations[i] - ground_truth[i];
    res = res.array() * res.array();
    rmse += res;
	}

	//calculating the mean
	rmse = rmse.array() / size;

	//calculating the squared root
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  MatrixXd Hj(3, 4);

	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

  float divOne = px * px + py * py;
	if (divOne == 0) {
    // create a simple matrix to manage first 2 vector values
    // (1, 0, 0, 0)
    // (0, 1, 0, 0)
    // (0, 0, 0, 0)
    Hj = MatrixXd::Identity(3, 4);
    Hj(2, 2) = 0;
    return Hj;
	}

  float divSqrt = sqrt(divOne);
  float divPow = divOne * divSqrt;

  float pxDividedBySqrt = px / divSqrt;
  float pyDividedBySqrt = py / divSqrt;

	Hj << pxDividedBySqrt, pyDividedBySqrt, 0, 0,
	    -py / divOne, px / divOne, 0, 0,
	    py * (vx * py - vy * px) / divPow, px * (vy * px - vx * py) / divPow, pxDividedBySqrt, pyDividedBySqrt;

	return Hj;
}

void Tools::CalculateCartesian(const float &rho, const float &phi, float &px, float &py) {
  px = rho * sin(phi);
  py = -rho * cos(phi);
}
