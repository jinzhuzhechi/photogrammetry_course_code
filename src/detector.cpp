#include "detector.h"

using namespace cv;
using namespace std;


vector<cv::Point> getRandomSubset(const vector<cv::Point>& input, int subsetSize) 
{
    std::vector<cv::Point> subset;
    // 创建随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, input.size() - 1);

    // 从输入向量中随机选择元素
    for (int i = 0; i < subsetSize; i++) {
        int randomIndex = dis(gen);
        subset.push_back(input[randomIndex]);
    }

    return subset;
}

vector<cv::Point> Harris_Corner_Detector(const cv::Mat &srcImg, double threshold, double distanceThreshold)
{
    vector<cv::Point> corners;
    vector<pair<float, cv::Point>> sortedCorners;

    cv::Mat imageGrey = srcImg.clone();
    if (imageGrey.channels() != 1)
        cv::cvtColor(imageGrey, imageGrey, cv::COLOR_BGR2GRAY);

    cv::Mat sobel_x = (cv::Mat_<int>(3, 3) << -1, 0, 1, -2, 0, 2, -1, 0, 1);
    cv::Mat sobel_y = (cv::Mat_<int>(3, 3) << -1, -2, -1, 0, 0, 0, 1, 2, 1);
    cv::Mat gauss_k = cv::getGaussianKernel(3, 0.5, CV_32F);

    cv::Mat dx, dy, dxdy;
    cv::filter2D(imageGrey, dx, -1, sobel_x);
    cv::filter2D(imageGrey, dy, -1, sobel_y);
    cv::multiply(dx, dy, dxdy);

    cv::Mat dx2, dy2;
    cv::multiply(dx, dx, dx2);
    cv::multiply(dy, dy, dy2);
    cv::filter2D(dx2, dx2, -1, gauss_k);
    cv::filter2D(dy2, dy2, -1, gauss_k);
    cv::filter2D(dxdy, dxdy, -1, gauss_k);

    double k = 0.04;

    cv::Mat harris = dx2.mul(dy2) - dxdy.mul(dxdy) - k * (dx2 + dy2).mul(dx2 + dy2);
    cv::normalize(harris, harris, 0, 255, cv::NORM_MINMAX, CV_32FC1);
    cv::convertScaleAbs(harris, harris);

    // 筛选响应值大于阈值的点作为角点
    for (int i = 100; i < harris.cols-100; i++)
        for (int j = 100; j < harris.rows-100; j++)
            if (harris.at<float>(i, j) > threshold)
                sortedCorners.emplace_back(harris.at<float>(i, j), cv::Point(i, j));

    
    // 按响应值排序
    sort(sortedCorners.begin(), sortedCorners.end(), [](const pair<float, cv::Point>& a, const pair<float, cv::Point>& b) {
        return a.first > b.first;
    });

    int n = std::min(50000, static_cast<int>(sortedCorners.size()));
    // 构造新的向量，包含前 n 个元素
    vector<pair<float, cv::Point>> subset(sortedCorners.begin(), sortedCorners.begin() + n);
    sortedCorners = subset;

    cout<<sortedCorners.size()<<endl;
    // 非极大值抑制
    vector<bool> survived(sortedCorners.size(), true);
    for (int i = 0; i < sortedCorners.size(); i++) {
        if (!survived[i]) continue;
        for (int j = i + 1; j < sortedCorners.size(); j++) {
            if (norm(sortedCorners[i].second - sortedCorners[j].second) < distanceThreshold) {
                survived[j] = false;
            }
        }
    }

    // 添加survived的角点
    for (int i = 0; i < sortedCorners.size(); i++) {
        if (survived[i]) {
            corners.push_back(sortedCorners[i].second);
        }
    }

    return corners;
}

vector<cv::Point> Moravec_Corner_Detector(const cv::Mat &srcImg, int windowSize, int threshold)
{
    vector<cv::Point> corners;

    if (windowSize % 2 == 0)
        windowSize += 1;

    cv::Mat img = srcImg.clone();
    if (img.channels() != 1)
        cv::cvtColor(img, img, cv::COLOR_RGB2GRAY);

    int r = windowSize / 2;
    cv::GaussianBlur(img, img, cv::Size(windowSize, windowSize), 0, 0);

    cv::Mat interestImg = cv::Mat::zeros(img.size(), CV_32FC1);

    /// 计算兴趣值
    for (int i = r; i < img.rows - r; i++)
        for (int j = r; j < img.cols - r; j++)
        {
            double value[4] = {0.0};
            double minValue = 0.0;

            for (int k = -r; k <= r; k++)
            {
                value[0] += pow(img.at<uchar>(i + k, j) - img.at<uchar>(i + k + 1, j), 2);
                value[1] += pow(img.at<uchar>(i, j + k) - img.at<uchar>(i, j + k + 1), 2);
                value[2] += pow(img.at<uchar>(i + k, j + k) - img.at<uchar>(i + k + 1, j + k + 1), 2);
                value[3] += pow(img.at<uchar>(i + k, j - k) - img.at<uchar>(i + k + 1, j - k - 1), 2);
            }

            minValue = min(min(value[0], value[1]), min(value[2], value[3]));
            interestImg.at<float>(i, j) = minValue;
        }

    /// 选取候选点
    int maxValue;
    cv::Point point;
    for (int i = r; i < img.rows - r;)
    {
        for (int j = r; j < img.cols - r;)
        {
            point.x = -1;
            point.y = -1;
            maxValue = 0;
            for (int m = -r; m < r; m++)
            {
                for (int n = -r; n < r; n++)
                {
                    if (interestImg.at<float>(i + m, j + n) > maxValue)
                    {
                        maxValue = interestImg.at<float>(i + m, j + n);
                        point.y = i + m;
                        point.x = j + n;
                    }
                }
            }
            if (maxValue > threshold)
            {
                corners.push_back(point);
            }
            j += windowSize;
        }
        i += windowSize;
    }

    return corners;
}

vector<cv::Point> SIFT_Detector(const cv::Mat &srcImg)
{
    vector<cv::Point> corners;

    vector<cv::KeyPoint> keypoints;
    Ptr<cv::FeatureDetector> detector = cv::SIFT::create();
    detector->detect(srcImg, keypoints);

    for (cv::KeyPoint kpt : keypoints)
        corners.push_back(cv::Point(kpt.pt));

    return corners;
}

void Draw_Corner(const cv::Mat &srcImg, cv::Mat &outputImg, vector<cv::Point> corners)
{
    outputImg = srcImg.clone();
    static default_random_engine e;
    static uniform_int_distribution<int> u(0, 255);

    for (cv::Point corner : corners)
    {
        cv::Scalar color(u(e), u(e), u(e));
        cv::circle(outputImg, corner, 5, color, 2);
    }
}