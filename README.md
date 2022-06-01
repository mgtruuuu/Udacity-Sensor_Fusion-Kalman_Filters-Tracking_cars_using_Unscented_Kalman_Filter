# Tracking cars using Unscented Kalman Filter

<img src="media/ukf_highway_tracked.gif" width="700" height="400" />

## 1. Project Overview

This project deals with **combining different sensor data with Kalman filters** to perceive and track objects over time. In this project Unscented Kalman Filter is used to estimate the state of multiple cars on a highway using noisy LiDAR and RaDAR measurements.


Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project rubric. 




<img src="media/ukf_highway.png" width="700" height="400" />


## Structure

`main.cpp` is using `highway.h` to create a straight 3 lane highway environment with 3 traffic cars and the main ego car at the center. The viewer scene is centered around the ego car and the coordinate system is relative to the ego car as well. The ego car is green while the other traffic cars are blue. The traffic cars will be accelerating and altering their steering to change lanes. Each of the traffic car's has it's own UKF object generated for it, and will update each individual one during every time step. 

The red spheres above cars represent the (x,y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.





## 2. Key Implementation (Unscented Kalman Filter)

![UKF-roadmap](/media/UKF-roadmap.PNG)



## 3. Dependencies for Running Locally


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




## 5. Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar modify the code in `highway.h` to alter the cars. Also check out `tools.cpp` to change how measurements are taken, for instance lidar markers could be the (x,y) center of bounding boxes by scanning the PCD environment and performing clustering. This is similar to what was done in the project [Sensor Fusion Lidar Obstacle Detection](https://github.com/mgtruuuu/Udacity-Sensor_Fusion-LiDAR-Lidar_Obstacle_Detection).
