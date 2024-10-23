#include "matcher.h"

using namespace std;
using namespace Eigen;
using namespace cv;

/// 设置匹配窗口大小
void BaseMatcher::setWindowSize(int windowSize)
{
    this->windowSize = windowSize;
}

/// 设置相关系数阈值
void BaseMatcher::setThreshold(double threshold)
{
    this->threshold = threshold;
}

/// 计算两个图像窗口的相关系数
double BaseMatcher::computeCorrelationIdx(const cv::Mat &srcWindow, const cv::Mat &dstWindow)
{
    if (srcWindow.size != dstWindow.size)
    {
        std::cerr << "窗口大小不匹配！" << std::endl;
        return 0;
    }
    // 总像元数量
    double totalNum = srcWindow.rows * srcWindow.cols;

    // 计算两张影像窗口中像素的平均值
    double gAverSrc = 0.0, gAverDst = 0.0;
    for (int i = 0; i < srcWindow.rows; i++)
        for (int j = 0; j < srcWindow.cols; j++)
        {
            gAverSrc += srcWindow.at<unsigned char>(i, j);
            gAverDst += dstWindow.at<unsigned char>(i, j);
        }
    gAverSrc /= totalNum;
    gAverDst /= totalNum;

    // 计算相关系数
    double numerator = 0.0;
    double denominatorSrc = 0.0;
    double denominatorDst = 0.0;

    for (int i = 0; i < srcWindow.rows; i++)
        for (int j = 0; j < srcWindow.cols; j++)
        {
            numerator += (srcWindow.at<unsigned char>(i, j) - gAverSrc) * (dstWindow.at<unsigned char>(i, j) - gAverDst);
            denominatorSrc += pow((srcWindow.at<unsigned char>(i, j) - gAverSrc), 2);
            denominatorDst += pow((dstWindow.at<unsigned char>(i, j) - gAverDst), 2);
        }

    double denominator = sqrt(denominatorSrc * denominatorDst);
    return abs(numerator / denominator);
}

/// 判断当前点是否在有效的图像范围之内
bool BaseMatcher::isVaildPoint(const cv::Mat &srcImg, const cv::Point &pt)
{
    double r = windowSize / 2;
    if (pt.x < r || pt.y < r || pt.x >= srcImg.cols - r || pt.y >= srcImg.rows - r)
        return false;
    return true;
}

/// 绘制同名点匹配结果
void BaseMatcher::drawMatches(const cv::Mat &srcImg, const cv::Mat &dstImg, cv::Mat &outputImg, std::vector<Match> &matches)
{
    outputImg.create(cv::Size(srcImg.cols + dstImg.cols, std::max(srcImg.rows, dstImg.rows)), CV_8UC3);

    srcImg.copyTo(outputImg(cv::Rect(0, 0, srcImg.cols, srcImg.rows)));
    dstImg.copyTo(outputImg(cv::Rect(srcImg.cols, 0, dstImg.cols, dstImg.rows)));

    cv::Point pt1, pt2;
    static std::default_random_engine e;
    static std::uniform_int_distribution<int> u(0, 255);

    for (Match match : matches)
    {
        cv::Scalar color(u(e), u(e), u(e));

        pt1 = match.srcPt;
        pt2 = cv::Point(match.dstPt.x + srcImg.cols, match.dstPt.y);

        cv::circle(outputImg, pt1, 5, color, 2);
        cv::circle(outputImg, pt2, 5, color, 2);
        cv::line(outputImg, pt1, pt2, color, 2);
    }
}



/// 同名点匹配，在不知道初始视差的情况下，使用相关系数匹配
void CorrelationMatcher::Correlation_Match(const cv::Mat &srcImg, 
                                            const cv::Mat &dstImg, 
                                            const std::vector<cv::Point> &srcPts, 
                                            const std::vector<cv::Point> &dstPts, 
                                            std::vector<Match> &matches)
{
    matches.clear();

    cv::Mat srcImgCopy, dstImgCopy;
    if (srcImg.channels() != 1)
        cv::cvtColor(srcImg, srcImgCopy, cv::COLOR_BGR2GRAY);
    else
        srcImgCopy = srcImg.clone();

    if (dstImg.channels() != 1)
        cv::cvtColor(dstImg, dstImgCopy, cv::COLOR_BGR2GRAY);
    else
        dstImgCopy = dstImg.clone();

    double r = windowSize / 2;
    int num = 0;

    int total = srcPts.size();

    for (cv::Point srcPt : srcPts)
    {
        if (isVaildPoint(srcImgCopy, srcPt))
        {
            // 使用ROI从原图中切出一个窗口
            cv::Rect rectSrc, rectDst;
            cv::Mat windowSrc, windowDst;
            rectSrc = cv::Rect(srcPt.x - r, srcPt.y - r, windowSize, windowSize);
            windowSrc = srcImgCopy(rectSrc);

            // 遍历目标图像的特征点，寻找满足条件的同名点
            double idx = 0.0, maxIdx = 0.0;
            cv::Point maxPt;
            for (cv::Point dstPt : dstPts)
            {
                if (isVaildPoint(dstImgCopy, dstPt))
                {
                    rectDst = cv::Rect(dstPt.x - r, dstPt.y - r, windowSize, windowSize);
                    windowDst = dstImgCopy(rectDst);
                    idx = computeCorrelationIdx(windowSrc, windowDst);
                    if (idx > maxIdx)
                    {
                        maxIdx = idx;
                        maxPt = dstPt;
                    }
                }
            }
            // 判断最大的相关系数是否满足设定阈值
            if (maxIdx > threshold)
            {
                Match match;
                match.srcPt = srcPt;
                match.dstPt = maxPt;
                match.dist = maxIdx;
                matches.push_back(match);
            }
        }
        num++;
    }
}
// 在已知初始视差的情况下，使用相关系数匹配
void CorrelationMatcher::Correlation_Match_initiated (const cv::Mat &srcImg_ori,
                                                        const cv::Mat &dstImg_ori,
                                                        const std::vector<cv::Point> &srcPts,
                                                        const int dx,
                                                        const int dy,
                                                        vector<Match> &matches)
{
    matches.clear();

    cv::Mat srcImg, dstImg;
    if (srcImg.channels() != 1)
        cv::cvtColor(srcImg_ori, srcImg, cv::COLOR_BGR2GRAY);
    else
        srcImg = srcImg_ori.clone();

    if (dstImg.channels() != 1)
        cv::cvtColor(dstImg_ori, dstImg, cv::COLOR_BGR2GRAY);
    else
        dstImg = dstImg_ori.clone();

    double r = windowSize / 2;
    int num = 0;

    int total = srcPts.size();

    for (cv::Point srcPt : srcPts)
    {
        if (isVaildPoint(srcImg, srcPt))
        {
            // 使用ROI从原图中切出一个窗口
            cv::Rect rectSrc, rectDst;
            cv::Mat windowSrc, windowDst;
            rectSrc = cv::Rect(srcPt.x - r, srcPt.y - r, windowSize, windowSize);
            windowSrc = srcImg(rectSrc);

            // 遍历目标图像的特征点，寻找满足条件的同名点
            double Coefficient = 0.0, maxCoefficient = 0.0;
            
            int x_start = max(srcPt.x + dx - 5,int(r)+1);
            int x_end = min(srcPt.x + dx + 5, dstImg.cols-int(r)-1);

            int y_start = max(srcPt.y + dy - 5,int(r)+1);
            int y_end = min(srcPt.y + dy + 5, dstImg.rows-int(r)-1);

            if (x_start >= x_end || y_start >= y_end)
                continue;

            int maxI, maxJ;
            for (int i = x_start; i < x_end; i++)
            {
                for (int j = y_start; j < y_end; j++)
                {
                    rectDst = cv::Rect(i - r, j - r, windowSize, windowSize);
                    windowDst = dstImg(rectDst);
                    Coefficient = computeCorrelationIdx(windowSrc, windowDst);
                    if (Coefficient > maxCoefficient)
                    {
                        maxCoefficient = Coefficient;
                        maxI = i;
                        maxJ = j;
                    }
                }
            }
            // 判断最大的相关系数是否满足设定阈值
            if (maxCoefficient > threshold)
            {
                Match match;
                match.srcPt = srcPt;
                match.dstPt = cv::Point2d(maxI, maxJ);
                match.dist = maxCoefficient;
                matches.push_back(match);
            }
        }
        num++;
    }
}

/// 单点最小二乘匹配
bool LsqMatcher::SubPixel_Match(const cv::Mat &srcImg, const cv::Mat &dstImg, Match &match)
{
    cv::Mat srcImgCopy, dstImgCopy;
    if (srcImg.channels() != 1)
        cv::cvtColor(srcImg, srcImgCopy, cv::COLOR_BGR2GRAY);
    else
        srcImgCopy = srcImg.clone();

    if (dstImg.channels() != 1)
        cv::cvtColor(dstImg, dstImgCopy, cv::COLOR_BGR2GRAY);
    else
        dstImgCopy = dstImg.clone();

    // 这里处理时认为x表示行，y表示列
    double y1 = match.srcPt.x;
    double x1 = match.srcPt.y;
    double y2 = match.dstPt.x;
    double x2 = match.dstPt.y;

    if (windowSize % 2 == 0)
        this->windowSize += 1;

    int r = windowSize / 2;

    // 使用ROI从原图中切出一个窗口
    cv::Rect rectSrc, rectDst;
    cv::Mat windowSrc, windowDst;
    rectSrc = cv::Rect(y1 - r, x1 - r, windowSize, windowSize);
    windowSrc = srcImgCopy(rectSrc);
    windowDst.create(cv::Size(windowSize, windowSize), CV_8UC1);

    // 设定几何畸变初值
    a0 = x2 - x1;
    a1 = 1;
    a2 = 0;
    b0 = y2 - y1;
    b1 = 0;
    b2 = 1;

    // 设定灰度畸变初值
    h0 = 0;
    h1 = 1;

    double xs = 0.0, ys = 0.0;
    double currentCorrelationIdx, bestCorrelationIdx = 0.0;
    cv::Point2d bestPt;

    for (int iter = 0; iter < 50; iter++) // 设定最大迭代次数不超过50次
    {
        Eigen::MatrixXd A(windowSize * windowSize, 8), L(windowSize * windowSize, 1), x;

        int num = 0;
        double xNumerator = 0.0, yNumerator = 0.0, xDenominator = 0.0, yDenominator = 0.0;

        for (int i = x1 - r; i <= x1 + r; i++)
            for (int j = y1 - r; j <= y1 + r; j++)
            {
                // 几何变形改正
                double m = a0 + a1 * i + a2 * j;
                double n = b0 + b1 * i + b2 * j;

                int I = floor(m);
                int J = floor(n);

                // 如果当前的点在图像的边界附近，就舍弃当前点，因为后面求导会出现问题
                if (I < 1 || I > dstImgCopy.rows || J < 1 || J > dstImgCopy.cols)
                    continue;

                // 重采样：双线性内插
                double pixelValue = (J + 1 - n) * ((I + 1 - m) * dstImgCopy.at<uchar>(I, J) + (m - I) * dstImgCopy.at<uchar>(I + 1, J)) + (n - J) * ((I + 1 - m) * dstImgCopy.at<uchar>(I, J + 1) + (m - I) * dstImgCopy.at<uchar>(I + 1, J + 1));

                // 辐射畸变改正
                pixelValue = h0 + h1 * pixelValue;
                windowDst.at<uchar>(i - x1 + r, j - y1 + r) = pixelValue;

                // 构建误差方程
                double gxDst = 0.5 * (dstImgCopy.at<uchar>(I + 1, J) - dstImgCopy.at<uchar>(I - 1, J));
                double gyDst = 0.5 * (dstImgCopy.at<uchar>(I, J + 1) - dstImgCopy.at<uchar>(I, J - 1));
                A(num, 0) = 1;
                A(num, 1) = pixelValue;
                A(num, 2) = gxDst;
                A(num, 3) = m * gxDst;
                A(num, 4) = n * gxDst;
                A(num, 5) = gyDst;
                A(num, 6) = m * gyDst;
                A(num, 7) = n * gyDst;

                L(num, 0) = srcImgCopy.at<uchar>(i, j) - pixelValue;

                // 计算最佳匹配点位
                double gxSrc = 0.5 * (srcImgCopy.at<uchar>(i + 1, j) - srcImgCopy.at<uchar>(i - 1, j));
                double gySrc = 0.5 * (srcImgCopy.at<uchar>(i, j + 1) - srcImgCopy.at<uchar>(i, j - 1));

                xNumerator += i * gxSrc * gxSrc;
                xDenominator += gxSrc * gxSrc;
                yNumerator += j * gySrc * gySrc;
                yDenominator += gySrc * gySrc;

                num++;
            }
        if (num < 8) // 无法求解法方程
            return false;

        currentCorrelationIdx = computeCorrelationIdx(windowSrc, windowDst);

        // std::cout << "Iter time: " << iter << std::endl;
        // std::cout << "a0 = " << a0 << "\ta1 = " << a1 << "\ta2 = " << a2 << std::endl;
        // std::cout << "b0 = " << b0 << "\tb1 = " << b1 << "\tb2 = " << b2 << std::endl;
        // std::cout << "h0 = " << h0 << "\th1 = " << h1 << std::endl;
        // std::cout << "idx = " << currentCorrelationIdx << std::endl;

        // std::cout << "A: \n" << A << std::endl;
        // std::cout << "L: \n" << L << std::endl;

        // 计算变形参数
        x = (A.transpose() * A).inverse() * (A.transpose() * L);
        // std::cout << "x: \n" << x << std::endl;
        // std::cout << std::endl;


        double a0_old = a0;
        double a1_old = a1;
        double a2_old = a2;
        double b0_old = b0;
        double b1_old = b1;
        double b2_old = b2;
        double h0_old = h0;
        double h1_old = h1;

        a0 = a0_old + x(2, 0) + a0_old * x(3, 0) + b0_old * x(4, 0);
        a1 = a1_old + a1_old * x(3, 0) + b1_old * x(4, 0);
        a2 = a2_old + a2_old * x(3, 0) + b2_old * x(4, 0);
        b0 = b0_old + x(5, 0) + a0_old * x(6, 0) + b0_old * x(7, 0);
        b1 = b1_old + a1_old * x(6, 0) + b1_old * x(7, 0);
        b2 = b2_old + a2_old * x(6, 0) + b2_old * x(7, 0);
        h0 = h0_old + x(0, 0) + h0_old * x(1, 0);
        h1 = h1_old + h1_old * x(1, 0);

        // 计算最佳匹配点位
        double xt = xNumerator / xDenominator;
        double yt = yNumerator / yDenominator;

        xs = a0 + a1 * xt + a2 * yt;
        ys = b0 + b1 * xt + b2 * yt;

        if (currentCorrelationIdx > bestCorrelationIdx)
        {
            bestPt.x = ys;
            bestPt.y = xs;
            bestCorrelationIdx = currentCorrelationIdx;
        }

        if (bestCorrelationIdx > threshold)
        {
            match.dstPt.x = bestPt.x;
            match.dstPt.y = bestPt.y;
            match.dist = bestCorrelationIdx;
            return true;
        }
    }

    match.dstPt.x = bestPt.x;
    match.dstPt.y = bestPt.y;
    match.dist = bestCorrelationIdx;
    return true;
}

cv::Mat Get_Homography_Matrix(const cv::Mat &srcImg_ori, const cv::Mat &dstImg_ori)
{
    cv::Mat srcImg, dstImg;
    if (srcImg_ori.channels() != 1)
        cv::cvtColor(srcImg_ori, srcImg, cv::COLOR_BGR2GRAY);
    else
        srcImg = srcImg_ori.clone();

    if (dstImg.channels() != 1)
        cv::cvtColor(dstImg_ori, dstImg, cv::COLOR_BGR2GRAY);
    else
        dstImg = dstImg_ori.clone();
    
    // 创建sift特征提取器

    Ptr<SIFT> sift = SIFT::create();

    vector<KeyPoint> keypoints1, keypoints2;
    Mat descriptors1, descriptors2;

    sift->detectAndCompute(srcImg, noArray(), keypoints1, descriptors1);
    sift->detectAndCompute(dstImg, noArray(), keypoints2, descriptors2);

    // 创建描述子匹配器
    Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create(DescriptorMatcher::FLANNBASED);
    vector<vector<DMatch>> knnMatches;
    matcher->knnMatch(descriptors1, descriptors2, knnMatches, 2);

    // 筛选匹配点
    vector<DMatch> goodMatches;
    float ratioThreshold = 0.7;  // 匹配距离的阈值

    for (size_t i = 0; i < knnMatches.size(); ++i) 
    {
        if (knnMatches[i][0].distance < ratioThreshold * knnMatches[i][1].distance) 
        {
            goodMatches.push_back(knnMatches[i][0]);
        }
    }
    // 获取匹配点的坐标
    vector<Point2f> points1, points2;

    for (const auto& match : goodMatches) 
    {
        points1.push_back(keypoints1[match.queryIdx].pt);
        points2.push_back(keypoints2[match.trainIdx].pt);
    }

    Mat homography = findHomography(points1, points2, RANSAC);
    return homography;
}