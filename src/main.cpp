/* \author Aaron Brown */
// Create simple 3d highway enviroment using PCL
// for exploring self-driving car sensors

//#include "render/render.h"
#include "highway.h"


std::vector<double> lidarNISs{};	// to check the distribution of NIS
std::vector<double> radarNISs{};	// to check the distribution of NIS


int main(int argc, char** argv) {

	pcl::visualization::PCLVisualizer::Ptr viewer{ new pcl::visualization::PCLVisualizer("3D Viewer") };
	viewer->setBackgroundColor(0, 0, 0);

	// Set camera position and angle.
	viewer->initCameraParameters();
	constexpr float x_pos{ .0f };
	viewer->setCameraPosition(x_pos - 26, 0, 15.0, x_pos + 25, 0, 0, 0, 0, 1);

	Highway highway{ viewer };

	//initHighway(viewer);

	constexpr int frame_per_sec{ 30 };
	constexpr int sec_interval{ 10 };
	int frame_count{ 0 };
	int time_us{ 0 };

	const double egoVelocity{ 25 };

	while (frame_count < (frame_per_sec * sec_interval)) {
		viewer->removeAllPointClouds();
		viewer->removeAllShapes();

		//stepHighway(egoVelocity,time_us, frame_per_sec, viewer);
		highway.stepHighway(egoVelocity, time_us, frame_per_sec, viewer);
		viewer->spinOnce(1000 / frame_per_sec);
		++frame_count;
		time_us = 1000000 * frame_count / frame_per_sec; 
	}

	
	// Exclude first 6 NIS values for all cars(3).
	int num_overLidarNIS = std::count_if(lidarNISs.begin() + (6 * 3), lidarNISs.end(), [](double e) { return e > 5.991; });
		
	// Exclude first 6 NIS values for all cars(3).
	int num_overRadarNIS = std::count_if(radarNISs.begin() + (6 * 3), radarNISs.end(), [](double e) { return e > 7.815; });

	std::cout << "The ratio of [# LiDAR NIS] over limit : " << static_cast<double>(num_overLidarNIS) / (lidarNISs.size() - (6 * 3));
	std::cout << "The ratio of [# RaDAR NIS] over limit : " << static_cast<double>(num_overRadarNIS) / (radarNISs.size() - (6 * 3));




	return 0;
}