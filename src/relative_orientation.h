#ifndef RELATIVE_ORIENTATION
#define RELATIVE_ORIENTATION


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <random>


#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"

#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"


#include "photo.h"

struct point_pair
{
    double x1, x2, y1, y2;
    point_pair(double x1_, double y1_, double x2_, double y2_) : x1(x1_), y1(y1_), x2(x2_), y2(y2_) {}
};

// 获取同名点对
vector<point_pair> get_corresponding_point_pair(cv::Mat &image1, cv::Mat &image2);

// 相对定向
void get_relative_orientation(cv::Mat &img1, cv::Mat &img2, double f);

#endif