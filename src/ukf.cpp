#include "ukf.h"

#include "Eigen/Dense"
//#include <Eigen/Dense>

#include <iostream>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

extern std::vector<double> lidarNISs;
extern std::vector<double> radarNISs;




// Given a yaw (angle) in rad, returns the mod value in the range (-pi, pi].
double mod_2Pi_MPi2Pi(double angle) {

	// Find x such that angle = x (mod 2 * Pi ) and -Pi < x <= +Pi.
	while (angle > M_PI)	angle -= 2. * M_PI;
	while (angle <= -M_PI)	angle += 2. * M_PI;

	return angle;
};




/**
 * Initializes Unscented Kalman filter.
 */
UKF::UKF() {

	// If this is false, laser measurements will be ignored (except during init).
	use_laser_ = true;

	// If this is false, radar measurements will be ignored (except during init).
	use_radar_ = true;

	// Initial state vector.
	x_ = VectorXd(5);

	//
	//// Initial covariance matrix.
	//
	P_ = MatrixXd::Identity(5, 5);
	//P_ = MatrixXd(5, 5);

	//
	//// Process noise standard deviation longitudinal acceleration in m/s^2
	//
	constexpr double max_acceleration{ 5.0 };
	std_a_ = max_acceleration / 2;

	//
	//// Process noise standard deviation yaw acceleration in rad/s^2
	//
	std_yawdd_ = 5.0;



   /**
	* DO NOT MODIFY measurement noise values below.
	* These are provided by the sensor manufacturer.
	*/

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
	* End
	* DO NOT MODIFY section for measurement noise values.
	*/





   /**
	* TODO: Complete the initialization. See ukf.h for other member properties.
	* Hint: one or more values initialized above might be wildly off...
	*/

	n_x_ = 5;
	n_aug_ = 7;
	n_sig_ = 2 * n_aug_ + 1;
	lambda_ = 3.0 - n_x_;

	Xsig_pred_ = MatrixXd(n_x_, n_sig_);
	Xsig_pred_.fill(0.0);
	




	////////????????????????? RIGHT NOW??

	// Set Weights to predict Mean and Covariance (Step 3)
	weights_ = VectorXd(n_sig_);
	weights_(0) = lambda_ / (lambda_ + n_aug_);
	for (int i{ 1 }; i < n_sig_; ++i)		// 2n+1 weights
		weights_(i) = 0.5 / (n_aug_ + lambda_);



	is_initialized_ = false;
}



UKF::~UKF() {}





void UKF::Prediction(double delta_t) {

	// STEP 1 : Generate Sigma Points.

	// Create augmented mean vector.
	VectorXd x_aug(n_aug_);
	x_aug.head(n_x_) = x_;				// Augment the mean state.
	x_aug(5) = 0;
	x_aug(6) = 0;

	// Create augmented state covariance.
	MatrixXd P_aug(n_aug_, n_aug_);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_x_, n_x_) = std_a_ * std_a_;
	P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

	const MatrixXd L{ P_aug.llt().matrixL() };	// Create square root matrix.

	// Create sigma point matrix.
	MatrixXd Xsig_aug(n_aug_, n_sig_);
	Xsig_aug.fill(0.0);
	Xsig_aug.col(0) = x_aug;			// Set first column of sigma point matrix.
	for (int i{ 0 }; i < n_aug_; ++i) {	// Set remaining sigma points.
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}

			

	// STEP 2 : Predict Sigma Points.

	const double dt{ delta_t };

	// Loop over each sigma point, transforming back to the original state space????
	for (int i{ 0 }; i < n_sig_; ++i) {

		// Extract values for better readability.
		const double p_x{ Xsig_aug(0, i) };
		const double p_y{ Xsig_aug(1, i) };
		const double v{ Xsig_aug(2, i) };
		const double yaw{ Xsig_aug(3, i) };
		const double yawd{ Xsig_aug(4, i) };
		const double nu_a{ Xsig_aug(5, i) };
		const double nu_yawdd{ Xsig_aug(6, i) };

		// predicted state values
		double v_p{ v };
		double yaw_p{ yaw + yawd * delta_t };
		double yawd_p{ yawd };
		double px_p, py_p;
		// Avoid division by zero.
		if (std::fabs(yawd) > 0.001) {
			px_p = p_x + (v / yawd) * (sin(yaw + yawd * dt) - sin(yaw));
			py_p = p_y + (v / yawd) * (-cos(yaw + yawd * dt) + cos(yaw));
		} else {
			px_p = p_x + v * dt * cos(yaw);
			py_p = p_y + v * dt * sin(yaw);
		}

		// Add noise.
		px_p += 0.5 * nu_a * delta_t * delta_t * cos(yaw);
		py_p += 0.5 * nu_a * delta_t * delta_t * sin(yaw);
		v_p += nu_a * delta_t;
		yaw_p += 0.5 * nu_yawdd * delta_t * delta_t;
		yawd_p += nu_yawdd * delta_t;

		// Write predicted sigma point into right column.
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}



	// STEP 3 : Predict Mean and Covariance.

	// Create vector for predicted state.
	VectorXd x(n_x_);
	x.fill(0.0);

	// Create covariance matrix for prediction.
	MatrixXd P(n_x_, n_x_);
	P.fill(0.0);

	// Predict the mean state.
	for (int i{ 0 }; i < n_sig_; ++i)
		x += (weights_(i) * Xsig_pred_.col(i));

	// Predict the covariance matrix.
	for (int i{ 0 }; i < n_sig_; ++i) {

		VectorXd x_diff{ Xsig_pred_.col(i) - x };		// state difference

		Xsig_pred_(3) = mod_2Pi_MPi2Pi(Xsig_pred_(3));	// angle normalization

		P += weights_(i) * x_diff * x_diff.transpose();
	}

	// Write result.
	x_ = x;
	P_ = P;
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

	// For the initial measurement, skip the prediction step.
	if (is_initialized_ == false) {

		// Use the first measurement to set the value of the mean state.
		x_.fill(0.0);
		x_.head(2) << meas_package.raw_measurements_;

		time_us_ = meas_package.timestamp_;

		is_initialized_ = true;
	}
	// For subsequent measurement, trigger a full prediction/update cycle.
	else {

		const double delta_t{ (meas_package.timestamp_ - time_us_) / 1e6 };

		time_us_ = meas_package.timestamp_;

		Prediction(delta_t);
	}

	if (MeasurementPackage::SensorType::LASER == meas_package.sensor_type_)
		UpdateLidar(meas_package);
	else if (MeasurementPackage::SensorType::RADAR == meas_package.sensor_type_)
		UpdateRadar(meas_package);
}



void UKF::UpdateLidar(MeasurementPackage meas_package) {

    // STEP 4 : case a) Predict LiDAR Measurements.

	const VectorXd z(meas_package.raw_measurements_);	// 2x1 vector for lidar
	const int n_z{ static_cast<int>(z.size())};

	// Create matrix in measurement space to make measurement sigma points.
	MatrixXd Zsig(n_z, n_sig_);
	Zsig = Xsig_pred_.topLeftCorner(n_z, n_sig_);

	// Create predicted measurement mean.
	VectorXd z_pred(n_z);
	z_pred.fill(0.0);
	for (int i{ 0 }; i < n_sig_; ++i)
		z_pred += weights_(i) * Zsig.col(i);

	// Create predicted measurement covariance.
	MatrixXd S(n_z, n_z);
	S.fill(0.0);
	for (int i{ 0 }; i < n_sig_; ++i) {
		VectorXd z_diff{ Zsig.col(i) - z_pred };	// residual
		z_diff(1) = mod_2Pi_MPi2Pi(z_diff(1));		// angle normalization

		S += weights_(i) * z_diff * z_diff.transpose();
	}
	MatrixXd R(n_z, n_z);	// Create measurement noise covariance.
	R << std_laspx_ * std_laspx_, 0,
		0, std_laspy_* std_laspy_;
	S += R;					// Add measurement noise covariance.
	


	// STEP 5) Update the state by applying the Kalman gain to the residual

	// Create matrix Tc for Cross-correlation 
	// between sigma points in state space and measurement space.
	MatrixXd Tc(n_x_, n_z);
	Tc.fill(0.0);
	for (int i{ 0 }; i < n_sig_; ++i) {
		VectorXd z_diff{ Zsig.col(i) - z_pred };	// residual
		z_diff(1) = mod_2Pi_MPi2Pi(z_diff(1));		// angle normalization

		VectorXd x_diff{ Xsig_pred_.col(i) - x_ };
		x_diff(3) = mod_2Pi_MPi2Pi(x_diff(3));		// angle normalization

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}

	// Update state mean and covariance.
	const MatrixXd K{ Tc * S.inverse() };			// Kalman gain K
	VectorXd residuals{ z - z_pred };				// residual
	residuals(1) = mod_2Pi_MPi2Pi(residuals(1));	// angle normalization
	x_ += K * residuals;
	P_ -= K * S * K.transpose();

	// Calculate normalized innovation squared (NIS) for tuning
	const double lidarNIS{ residuals.transpose() * S.inverse() * residuals };
	//cout << "\tLidar NIS ( ~ X^2 : P(e<5.991) = 0.95 for 2DF ) : " << lidarNIS << endl;
	lidarNISs.push_back(lidarNIS);
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

	// STEP 4 : case b) Predict RaDAR Measurements.

	const VectorXd z{ meas_package.raw_measurements_ };	// 3x1 vector for radar
	const int n_z{ static_cast<int>(z.size()) };			

	// Create matrix in measurement space to make measurement sigma points.
	MatrixXd Zsig(n_z, n_sig_);

	// Create predicted measurement mean.
	VectorXd z_pred(n_z);
	
	// Create predicted measurement covariance.
	MatrixXd S(n_z, n_z);

	// Transform sigma points into measurement space : predictedSP -> measurementSP
	for (int i{ 0 }; i < n_sig_; ++i) {

		const double p_x{ Xsig_pred_(0, i) };
		const double p_y{ Xsig_pred_(1, i) };
		const double v{ Xsig_pred_(2, i) };	
		const double yaw{ Xsig_pred_(3, i) };
		const double v1{ cos(yaw) * v };
		const double v2{ sin(yaw) * v };

		// Apply measurement model.
		Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);		// r
		Zsig(1, i) = atan2(p_y, p_x);					// phi
		Zsig(2, i) = (Zsig(0, i) > 0.001) ? (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y) : 0;	// r_dot
	}


	// Calculate predicted measurement mean.
	z_pred.fill(0.0);
	for (int i{ 0 }; i < n_sig_; ++i)
		z_pred += weights_(i) * Zsig.col(i);


	// Calculate predicted measurement covariance.
	S.fill(0.0);
	for (int i{ 0 }; i < n_sig_; ++i) {
		VectorXd z_diff{ Zsig.col(i) - z_pred };	// residual
		z_diff(1) = mod_2Pi_MPi2Pi(z_diff(1));		// angle normalization

		S += weights_(i) * z_diff * z_diff.transpose();
	}
	MatrixXd R(n_z, n_z);	// Create measurement noise covariance.
	R << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_* std_radphi_, 0,
		0, 0, std_radrd_* std_radrd_;
	S += R;					// Add measurement noise covariance.



	// STEP 5 : Update the state by applying the Kalman gain to the residual.

	// Create matrix Tc for Cross-correlation 
	// between sigma points in state space and measurement space.
	MatrixXd Tc(n_x_, n_z);
	Tc.fill(0.0);
	for (int i{ 0 }; i < n_sig_; ++i) {
		VectorXd z_diff{ Zsig.col(i) - z_pred };	// residual
		z_diff(1) = mod_2Pi_MPi2Pi(z_diff(1));		// angle normalization

		VectorXd x_diff{ Xsig_pred_.col(i) - x_ };
		x_diff(3) = mod_2Pi_MPi2Pi(x_diff(3));		// angle normalization

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}

	// Update state mean and covariance.
	const MatrixXd K{ Tc * S.inverse() };			// Kalman gain K
	VectorXd residuals{ z - z_pred };				// residual
	residuals(1) = mod_2Pi_MPi2Pi(residuals(1));	// angle normalization
	x_ += K * residuals;					
	P_ -= K * S * K.transpose();			

	// Calculate normalized innovation squared (NIS) for tuning
	const double radarNIS{ residuals.transpose() * S.inverse() * residuals };
	//cout << "Radar NIS ( ~ X^2 : P(e<7.815) = 0.95 for 3DF ) : " << radarNIS << endl;
	radarNISs.push_back(radarNIS);
}