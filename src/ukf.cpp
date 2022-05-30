#include "ukf.h"

#include "Eigen/Dense"
//#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;





// Given a yaw angle (psi) in rad, returns the mod value in the range (-pi, pi].
double mod_2Pi_MPi2Pi(double angle) {

	// Find x such that angle = x (mod 2 * Pi ) and -Pi < x <= +Pi.
	while (angle > M_PI)	angle -= 2. * M_PI;
	while (angle <= -M_PI)	angle += 2. * M_PI;

	return angle;
};




/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 30;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 30;

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
	* End DO NOT MODIFY section for measurement noise values
	*/

   /**
	* TODO: Complete the initialization. See ukf.h for other member properties.
	* Hint: one or more values initialized above might be wildly off...
	*/

	n_x_ = 5;
	n_aug_ = 7;
	n_sig_ = 2 * n_aug_ + 1;
	weights_ = VectorXd(n_sig_);
	lambda_ = 3.0 - n_x_;

	Xsig_pred_ = MatrixXd(n_x_, n_sig_);
	Xsig_pred_.fill(0.0);

	// Setup augmented weights vector
	double w0{ lambda_ / (lambda_ + n_aug_) };
	double w{ 1 / (2 * (lambda_ + n_aug_)) };
	weights_.fill(w);
	weights_(0) = w0;

	is_initialized_ = false;
}

UKF::~UKF() {}


// For each UKF instance initialize it.
void UKF::InitializeUKF(MeasurementPackage meas_package) {

	// Use the first measurement to set the value of the mean state.
	x_.fill(0.0);
	x_.head(2) << meas_package.raw_measurements_;

	time_us_ = meas_package.timestamp_;

	is_initialized_ = true;
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
   /**
	* TODO: Complete this function! Make sure you switch between lidar and radar
	* measurements.
	*/

	if (!is_initialized_) {

		// For the initial measurement, skip the prediction step
		InitializeUKF(meas_package);

		if (MeasurementPackage::SensorType::LASER == meas_package.sensor_type_)
			UpdateLidar(meas_package);
		else if (MeasurementPackage::SensorType::RADAR == meas_package.sensor_type_)
			UpdateRadar(meas_package);
	}
	else {

		// On each subsequent measurement, trigger a full prediction/update cycle
		double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;
		time_us_ = meas_package.timestamp_;
		Prediction(delta_t);

		if (MeasurementPackage::SensorType::LASER == meas_package.sensor_type_)
			UpdateLidar(meas_package);
		else if (MeasurementPackage::SensorType::RADAR == meas_package.sensor_type_)
			UpdateRadar(meas_package);
	}
}

void UKF::Prediction(double delta_t) {
   /**
	* TODO: Complete this function! Estimate the object's location.
	* Modify the state vector, x_. Predict sigma points, the state,
	* and the state covariance matrix.
	*/

	// STEP 1) Generate sigma (sampling) points in the augmented state space
// cout << "Generating sigma points in the augmented state space" << endl;

	MatrixXd Xsig = MatrixXd(n_aug_, n_sig_);
	Xsig.fill(0.0);
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;  // Augment the mean state

	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);  // Augment the covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_x_, n_x_) = pow(std_a_, 2);
	P_aug(n_x_ + 1, n_x_ + 1) = pow(std_yawdd_, 2);

	MatrixXd A_aug = P_aug.llt().matrixL();  // Create the square root matrix A
	double c_aug = sqrt(lambda_ + n_aug_);
	MatrixXd cA_aug = c_aug * A_aug;

	Xsig.col(0) = x_aug; // First column of Xsig is central sigma point
	for (int i = 1; i <= n_aug_; i++) {
		Xsig.col(i) = x_aug + cA_aug.col(i - 1);  // First group of sigma points
	}
	for (int i = n_aug_ + 1; i <= n_sig_ - 1; i++) {
		Xsig.col(i) = x_aug - cA_aug.col(i - 1 - n_aug_); // Symmetric group
	}

	// STEP 2) Predict the motion of each sigma point
	// cout << "Predicting the motion of each sigma point" << endl; 

	double dt = delta_t;
	Xsig_pred_.fill(0.0);

	// Loop over each sigma point, transforming back to the original state space
	for (int i = 0; i < n_sig_; i++) {
		double px = Xsig(0, i);
		double py = Xsig(1, i);
		double v = Xsig(2, i);
		double psi = Xsig(3, i);
		double psid = Xsig(4, i);
		double nu_a = Xsig(5, i);
		double nu_psidd = Xsig(6, i);

		// Transform from the augmented to the original (n_x-dim) state space

		// Avoid division by zero
		if (std::fabs(psid) > 0.001) {
			Xsig_pred_(0, i) = px + v / psid * (sin(psi + psid * dt) - sin(psi)) +
				pow(dt, 2) / 2 * cos(psi) * nu_a;
			Xsig_pred_(1, i) = py + v / psid * (-cos(psi + psid * dt) + cos(psi)) +
				pow(dt, 2) / 2 * sin(psi) * nu_a;
		}
		else {
			Xsig_pred_(0, i) = px + v * dt * cos(psi) +
				pow(dt, 2) / 2 * cos(psi) * nu_a;
			Xsig_pred_(1, i) = py + v * dt * sin(psi) +
				pow(dt, 2) / 2 * sin(psi) * nu_a;
		}

		Xsig_pred_(2, i) = v + 0 + dt * nu_a;
		Xsig_pred_(3, i) = psi + psid * dt + pow(dt, 2) / 2 * nu_psidd;
		Xsig_pred_(4, i) = psid + 0 + dt * nu_psidd;
	}

	// STEP 3) Predict the next state: mean and covariance matrices
	// cout << "Predicting the next state's mean and covariance" << endl;

	VectorXd x = VectorXd(n_x_);
	x.fill(0.0);
	MatrixXd P = MatrixXd(n_x_, n_x_);
	P.fill(0.0);

	// Predict the mean state
	for (int i = 0; i < n_sig_; i++) {
		x = x + weights_(i) * Xsig_pred_.col(i);
	}

	// Predict the covariance matrix
	for (int i = 0; i < n_sig_; i++) {
		P = P + weights_(i) *
			(Xsig_pred_.col(i) - x_) *
			(Xsig_pred_.col(i) - x_).transpose();
	}

	x_ = x;
	P_ = P;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
   /**
	* TODO: Complete this function! Use lidar data to update the belief
	* about the object's position. Modify the state vector, x_, and
	* covariance, P_.
    * You can also calculate the lidar NIS, if desired.
    */

	// STEP 4) Predict the measurement mean and covariance; calculate Kalman gain
// cout << "Predicting the measurement mean and covariance; calculating Kalman gain" << endl; 

	const VectorXd z{ meas_package.raw_measurements_ };
	const int n_z{ static_cast<int>(z.size())};	// Measurement z is a 2x1 vector for lidar

	// Measurement matrix
	MatrixXd H{ n_z, n_x_ };
	H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;

	// Measurement covariance matrix
	MatrixXd R{ n_z, n_z };
	R << std_laspx_ * std_laspx_, 0,
		 0, std_laspy_* std_laspy_;

	VectorXd z_pred(n_z);
	z_pred = x_.head(n_z); // Extract the px, py values from the state

	const VectorXd y{ z - z_pred };  // Calculate the residuals vector y
	const MatrixXd H_t{ H.transpose() };
	const MatrixXd S{ H * P_ * H_t + R };
	const MatrixXd K{ P_ * H_t * S.inverse() };

	// STEP 5) Update the state, by applying the Kalman gain to the residual
	// cout << "Updating the state by applying the Kalman gain to the residual" << endl; 

	MatrixXd I{ MatrixXd::Identity(n_x_, n_x_) };

	// Update the mean and covariance matrix
	x_ += (K * y);
	P_ = (I - K * H) * P_;

	// Calculate normalized innovation squared (NIS) for tuning
	// double NIS = y.transpose() * Sinv * y;
	// cout << "Lidar NIS (2-df X^2, 95% < 5.991) = " << NIS << endl; 
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

   /**
	* TODO: Complete this function! Use radar data to update the belief
	* about the object's position. Modify the state vector, x_, and
	* covariance, P_.
	* You can also calculate the radar NIS, if desired.
	*/

	// STEP 4) Predict the measurement mean (z_pred) and innovation covariance (S)
	// cout << "Predicting the measurement mean and covariance; calculating Kalman gain" << endl; 

	const VectorXd z{ meas_package.raw_measurements_ };
	const int n_z{ static_cast<int>(z.size()) };			// Measurement z is a 3x1 vector for radar

	MatrixXd Zsig{ n_z, n_sig_ };
	VectorXd z_pred(n_z);
	z_pred.fill(0.0);

	// Transform the sigma points into measurement space
	for (int i{ 0 }; i < n_sig_; ++i) {
		const double px{ Xsig_pred_(0, i) };
		const double py{ Xsig_pred_(1, i) };
		const double v{ Xsig_pred_(2, i) };
		const double psi{ Xsig_pred_(3, i) };
		const double psid{ Xsig_pred_(4, i) };

		const double rho{ sqrt(pow(px, 2) + pow(py, 2)) };
		const double phi{ std::atan2(py, px) };
		double rhod{ 0.0 };
		if (std::fabs(rho) > 0.001)
			rhod = (px * cos(psi) * v + py * sin(psi) * v) / rho;

		Zsig(0, i) = rho;
		Zsig(1, i) = phi;
		Zsig(2, i) = rhod;
	}

	// Calculate the predicted mean measurement z_pred
	for (int i{ 0 }; i < n_sig_; ++i)
		z_pred += weights_(i) * Zsig.col(i);


	// Calculate the innovation covariance matrix S
	MatrixXd S{ n_z, n_z };


	//// Given a yaw angle (psi) in rad, returns the mod value in the range (-pi, pi].
	//auto mod_2Pi_MPi2Pi{
	//	[](double angle) {
	//		// Find x such that angle = x (mod 2 * Pi ) and -Pi < x <= +Pi.
	//		while (angle > M_PI)	angle -= 2. * M_PI;
	//		while (angle <= -M_PI)	angle += 2. * M_PI;

	//		return angle;
	//	}
	//};

	S.fill(0.0);
	for (int i{ 0 }; i < n_sig_; ++i) {
		VectorXd z_diff{ Zsig.col(i) - z_pred };
		z_diff(1) = mod_2Pi_MPi2Pi(z_diff(1));

		S += weights_(i) * z_diff * z_diff.transpose();
	}

	//double mod_angle = 0.0;	?????????????

	// Add measurement noise covariance matrix.
	MatrixXd R{ n_z, n_z };
	R << std_radr_ * std_radr_, 0, 0,
		0, std_radphi_* std_radphi_, 0,
		0, 0, std_radrd_* std_radrd_;
	S += R;




	// STEP 5) Update the state, by applying the Kalman gain to the residual
	// cout << "Updating the state by applying the Kalman gain to the residual" << endl; 

	// Calculate the cross correlation matrix Tc
	MatrixXd Tc{ n_x_, n_z };
	Tc.fill(0.0);

	for (int i{ 0 }; i < n_sig_; ++i) {
		VectorXd z_diff{ Zsig.col(i) - z_pred };
		z_diff(1) = mod_2Pi_MPi2Pi(z_diff(1));

		VectorXd x_diff{ Xsig_pred_.col(i) - x_ };
		x_diff(3) = mod_2Pi_MPi2Pi(x_diff(3));

		Tc += weights_(i) * x_diff * z_diff.transpose();
	}


	
	// Update the mean and covariance matrix.
	MatrixXd K{ Tc * S.inverse() };			// Kalman gain
	VectorXd residuals{ z - z_pred };		// residual
	residuals(1) = mod_2Pi_MPi2Pi(residuals(1));
	x_ += K * residuals;
	P_ -= K * S * K.transpose();

	// Calculate normalized innovation squared (NIS) for tuning
	// double NIS = residuals.transpose() * Sinv * residuals;
	// cout << "Radar NIS (3-df X^2, 95% < 7.815) = " << NIS << endl;
}