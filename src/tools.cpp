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
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
	VectorXd diff_square(4);
	diff_square << 0, 0, 0, 0;
	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i) {
		// ... your code here
		VectorXd diff = estimations[i] - ground_truth[i];
		VectorXd diff_square_i = diff.array()*diff.array();
		diff_square += diff_square_i;
	}

	//calculate the mean
	// ... your code here
	VectorXd diff_s_norm = diff_square / estimations.size();

	//calculate the squared root
	// ... your code here
	rmse = diff_s_norm.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj_ = MatrixXd(3, 4);

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (c1==0) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
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
