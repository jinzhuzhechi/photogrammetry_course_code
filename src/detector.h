#ifndef DETECTOR
#define DETECTOR


#include <vector>
#include <random>

#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"

using namespace cv;
using namespace std;

vector<cv::Point> Harris_Corner_Detector(const cv::Mat &srcImg, double threshold = 120, double distanceThreshold = 70);
vector<cv::Point> Moravec_Corner_Detector(const cv::Mat &srcImg, int windowSize = 5, int threshold = 700);
vector<cv::Point> SIFT_Detector(const cv::Mat &srcImg);
void Draw_Corner(const cv::Mat &srcImg, cv::Mat &outputImg, std::vector<cv::Point> corners);

#endif