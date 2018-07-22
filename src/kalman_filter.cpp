#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.14159265

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the mean states and state covariance matrix
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_ * x_; //estimated prediction
	VectorXd y = z - z_pred; //error vector
	MatrixXd Ht = H_.transpose(); 
	MatrixXd S = H_ * P_ * Ht + R_; 
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si; // kalman gain

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);


	VectorXd z_pred = VectorXd(3);

	float p_pred = sqrt(px*px + py * py); //predicted distance
	float phi_pred = atan2(py, px); // predicted polar angle
	// predict velocity
	float v_pred; 
	if (p_pred == 0) {
		//check division by zero
		cout << " division by zero" << endl;
		return;
	}
	else {
		v_pred = (px * vx + py * vy) / p_pred;
	}

	//update preidiction vector
	z_pred << p_pred, phi_pred, v_pred;

	//error vector
	VectorXd y = z - z_pred;

	// check that polar angle is between +-pi
	while (y[1] > PI) { 
		y[1] -= 2 * PI;
	}
	while (y[1]  < -PI) {
		y[1] += 2 * PI;
	}

	
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si; //kalman gain

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
