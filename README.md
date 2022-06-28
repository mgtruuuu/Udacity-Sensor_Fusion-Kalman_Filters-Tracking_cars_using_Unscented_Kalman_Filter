# Tracking cars using Unscented Kalman Filter



## 1. Project Overview

![ukf_highway](media/ukf_highway.png)

This project deals with **combining different sensor data with Kalman filters** to perceive and track objects over time. In this project **Unscented Kalman Filter** is used to estimate the state of multiple cars on a highway using noisy LiDAR and RaDAR measurements and the **CTRV motion model** is used, which assumes constant velocity and turning rate.

The accuracy is being calculated by the **RMSE**, Root Mean Squared Error over each time step and for each car. The RMSE is looking at the error difference between the ground truth and estimated position and velocities. All the cars are considered in the RMSE value so if one car is particularly inaccurate it would affect the overall RMSE in this case. In the animation above the 4 RMSE values are X, Y, Vx, Vy from top to bottom.

![main](media/ukf_highway_tracked_case1.gif)

`main.cpp` is using `highway.h` to create a straight 3 lane highway environment with 3 traffic cars and the main ego car at the center. The viewer scene is centered around the ego car and the coordinate system is relative to the ego car as well. The ego car is green while the other traffic cars are blue. The traffic cars will be accelerating and altering their steering to change lanes. Each of the traffic car's has it's own UKF object generated for it, and will update each individual one during every time step. 

The red spheres above cars represent the (x,y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.


### Additional Project Details

In `highway.h` there are a number of parameters that can be modified to help with testing and understanding.

```c++
// Parameters 
// --------------------------------
// Set which cars to track with UKF
const std::vector<bool> trackCars{ true, true, true };
// Visualize sensor measurements
const bool visualize_lidar{ true };
const bool visualize_radar{ true };
const bool visualize_pcd{ false };
// Predict path in the future using UKF
const double projectedTime{ 0 };
const int projectedSteps{ 0 };
// --------------------------------
```

- The `trackCars` list can be used to toggle on/off cars for the UKF objects to track. The default is to track all three cars on the road, but for testing it can be nice to toggle to only track one at a time. For instance to only track the first car `{ true, false, false }`.

    ![Projecting the UKF path 2 seconds into the future](/media/ukf_highway_tracked_case2.gif)

- The animation above shows what it looks like if the `projectedTime` and `projectedSteps` variables are used. Also in the animation above the visualization for lidar and radar has been set to `false` and `projectedTime = 2` and `projectedSteps = 6`. The green spheres then show the predicted position for the car in the future over a 2 second interval. The number of steps increases the number of positions to interpolate the time interval.

    ![Visualizing high resolution point cloud data for traffic](/media/ukf_highway_tracked_case3.gif)

- The `visualize_pcd` parameter can be set to true to visualize the cars from the lidar's perspective from point clouds. The traffic pcd data is available in `src/sensors/data/pcd` directory. As a bonus assignment this data could be used to cluster the individual cars using techniques from Lidar Obstacle Detection and then bounding boxes could be fitted around the car clusters. The bounding box center (x,y) point could then be used to represent the lidar marker that would be fed into the UKF instead of the project pre-generated lidar marker from `tools.cpp` `lidarSense` function (This is similar to what was done in the project [Sensor Fusion Lidar Obstacle Detection](https://github.com/mgtruuuu/Udacity-Sensor_Fusion-LiDAR-Lidar_Obstacle_Detection)).





## 2. Key Implementation (Unscented Kalman Filter)

![UKF-roadmap](/media/UKF-roadmap.PNG)

### a) Prediction

```c++
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
```

### b) Update

```c++
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

	// Create matrix Tc for Cross-correlation between sigma points in state space and measurement space.
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
	cout << "Radar NIS ( ~ X^2 : P(e<7.815) = 0.95 for 3DF ) : " << radarNIS << endl;
}
```








## 3. Dependencies for Running Locally

NOTE : As of now this project has an error when running the .exe file on Windows environment. I've also checked that it's working properly on Linux. I'm still trying to figure out what the problem is. I'll fix the error soon.

1. PCL 1.2 : refer to the [official instructions](https://pointclouds.org/downloads/)

    - linux : $ `sudo apt install libpcl-dev`

    - windows : PS> `.\vcpkg.exe install pcl[visualization]:x64-windows`

    - mac : $ `brew install pcl`


2. cmake >= 3.5

    - All OSes: [click here for installation instructions](https://cmake.org/install/)


3. make >= 4.1 (Linux, Mac), 3.81 (Windows)

    - Linux: make is installed by default on most Linux distros

    - Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)

    - Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)


4. gcc/g++ >= 5.4

    - Linux: gcc / g++ is installed by default on most Linux distros

    - Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)

    - Windows: recommend using [MinGW](http://www.mingw.org/)




## 4. Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./ukf_highway`
