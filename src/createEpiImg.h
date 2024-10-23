#ifndef CREATEEPIIMG
#define CREATEEPIIMG

#include <iostream>
#include <string>
#include <vector>
#include <fstream>


#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"

#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"


#include "photo.h"

using namespace std;
using namespace cv;

void calculate_epiImg_EOP(EOP left_eop, EOP right_eop, EOP& left_epiImg_eop, EOP& right_epiImg_eop);
void calculate_M(EOP eop, EOP epiImg_eop, cv::Mat& M);
void calculate_epiImgSize(IOP iop,cv::Mat left_M,  int& epi_width, int& epi_height);
void creat_epiImg(cv::Mat img, cv::Mat & epiImg, IOP iop, cv::Mat M);
void crop_epiImg(cv::Mat& epiImg);





#endif
