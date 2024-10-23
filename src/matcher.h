#ifndef MATCHER
#define MATCHER


#include <random>
#include <vector>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"

#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"


using namespace std;
using namespace Eigen;
using namespace cv;

struct Match
{
    cv::Point2d srcPt; // 左片像点坐标
    cv::Point2d dstPt; // 右片像点坐标
    double dist;       // 相似性测度的计算值
};

class BaseMatcher
{
public:
    void setWindowSize(int windowSize);
    void setThreshold(double threshold);
    void drawMatches(const cv::Mat &srcImg,
                        const cv::Mat &dstImg,
                        cv::Mat &outputImg,
                        vector<Match> &matches);

protected:
    int windowSize = 15;    // 窗口大小
    double threshold = 0.8; // 阈值
    double computeCorrelationIdx(const cv::Mat &srcWindow, const cv::Mat &dstWindow);
    bool isVaildPoint(const cv::Mat &srcImg, const cv::Point &pt);
};

class CorrelationMatcher : public BaseMatcher
{
public:
    void Correlation_Match (const cv::Mat &srcImg,
                const cv::Mat &dstImg,
                const std::vector<cv::Point> &srcPts,
                const std::vector<cv::Point> &dstPts,
                vector<Match> &matches);
    void Correlation_Match_initiated (const cv::Mat &srcImg,
                const cv::Mat &dstImg,
                const std::vector<cv::Point> &srcPts,
                const int dx,
                const int dy,
                vector<Match> &matches);
};

class LsqMatcher : public BaseMatcher
{
public:
    bool SubPixel_Match(const cv::Mat &srcImg, const cv::Mat &dstImg, Match &match);

private:
    double a0 = 0, a1 = 1, a2 = 0;
    double b0 = 0, b1 = 0, b2 = 0;
    double h0 = 0, h1 = 1;
};

cv::Mat Get_Homography_Matrix(const cv::Mat &srcImg, const cv::Mat &dstImg);


#endif