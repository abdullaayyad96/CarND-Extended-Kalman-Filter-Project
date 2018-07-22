#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
	//initialize vector
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	
	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() == 0 || estimations.size() != ground_truth.size()) {
		cout << "Error with estimated states size for RMSE calculation" << endl;
		return rmse;
	}
	
	//accumulate squared residuals
	VectorXd diff_square(4);
	diff_square << 0, 0, 0, 0;
	for (int i = 0; i < estimations.size(); ++i) {
		VectorXd diff = estimations[i] - ground_truth[i];
		VectorXd diff_square_i = diff.array()*diff.array();
		diff_square += diff_square_i;
	}

	//calculate the mean
	VectorXd diff_s_norm = diff_square / estimations.size();

	//calculate the squared root
	rmse = diff_s_norm.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
	MatrixXd Hj_ = MatrixXd::Zero(3, 4);

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute some terms to avoid repeated calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (c1==0) {
		cout << "Error of division by Zero for jacobian calculation" << endl;
		return Hj_;
	}
	else {
		//compute the Jacobian matrix
		Hj_ << (px / c2), (py / c2), 0, 0,
			-(py / c1), (px / c1), 0, 0,
			py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;

	}

	return Hj_;
}
