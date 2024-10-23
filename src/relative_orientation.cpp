#include "relative_orientation.h"



vector<point_pair> get_corresponding_point_pair(cv::Mat &image1, cv::Mat &image2)
{
    // 初始化 SIFT 特征检测器
    cv::Ptr<cv::SIFT> sift = cv::SIFT::create();
    // 检测并计算特征点
    std::vector<cv::KeyPoint> keypoints1, keypoints2;
    cv::Mat descriptors1, descriptors2;
    sift->detectAndCompute(image1, cv::Mat(), keypoints1, descriptors1);
    sift->detectAndCompute(image2, cv::Mat(), keypoints2, descriptors2);

    // 使用 FLANN 匹配器进行特征点匹配
    cv::Ptr<cv::DescriptorMatcher> matcher = cv::DescriptorMatcher::create(cv::DescriptorMatcher::FLANNBASED);
    std::vector<std::vector<cv::DMatch>> matches;
    matcher->knnMatch(descriptors1, descriptors2, matches, 2);

    // 对匹配结果进行筛选,使用 Lowe 的 Ratio Test
    std::vector<cv::DMatch> goodMatches;
    for (size_t i = 0; i < matches.size(); i++)
    {
        if (matches[i][0].distance < 0.7 * matches[i][1].distance)
        {
            goodMatches.push_back(matches[i][0]);
        }
    }

    // 获取同名点的坐标
    std::vector<cv::Point2f> points1, points2;
    for (size_t i = 0; i < goodMatches.size(); i++)
    {
        points1.push_back(keypoints1[goodMatches[i].queryIdx].pt);
        points2.push_back(keypoints2[goodMatches[i].trainIdx].pt);
    }

    // 使用 RANSAC 算法进行错误点剔除
    std::vector<uchar> inliers(points1.size(), 0);
    cv::Mat fundamental = cv::findFundamentalMat(points1, points2, inliers, cv::RANSAC);

    // 更新同名点的坐标
    std::vector<cv::Point2f> inlierPoints1, inlierPoints2;
    for (size_t i = 0; i < inliers.size(); i++)
    {
        if (inliers[i])
        {
            inlierPoints1.push_back(points1[i]);
            inlierPoints2.push_back(points2[i]);
        }
    }

    // 获取同名点的坐标
    vector<point_pair> a_point_pair;
    for (size_t i = 0; i < inlierPoints1.size(); i++)
    {
        a_point_pair.push_back(point_pair(inlierPoints1[i].x, inlierPoints1[i].y, inlierPoints2[i].x, inlierPoints2[i].y));
    }

    return a_point_pair;
}

void get_relative_orientation(cv::Mat &img1, cv::Mat &img2, double f)
{
    vector<point_pair> match = get_corresponding_point_pair(img1, img2);
    int size = match.size();
    if (match.size() > 100)
    {
        size = 100;
        // 创建一个随机数生成器
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, static_cast<int>(match.size() - 1));

        // 随机选择 100 个元素
        vector<point_pair> new_vector;
        while (new_vector.size() < 100) {
            int index = dis(gen);
            new_vector.push_back(match[index]);
            match.erase(match.begin() + index);
        }
        // 现在 new_vector 包含了原 vector 的 100 个随机元素
        match = std::move(new_vector);
    }

    Eigen::VectorXd delta(5);
    int n = 0;
    double u, v, phi, omega, kappa;
    u = v = phi = omega = kappa = 0;
    Eigen::Vector3d X1, X2;
    Eigen::MatrixXd R(3, 3);
    Eigen::MatrixXd A(size, 5), L(size, 1);
    double Bx = match[0].x1 - match[0].x2;
    do {
        double By = Bx * u;
        double Bz = Bx * v;
        R(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
        R(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
        R(0, 2) = -sin(phi) * cos(omega);
        R(1, 0) = cos(omega) * sin(kappa);
        R(1, 1) = cos(omega) * cos(kappa);
        R(1, 2) = -sin(omega);
        R(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
        R(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
        R(2, 2) = cos(phi) * cos(omega);
        for (int i = 0; i < size; i++)
        {
            X1 = Eigen::Vector3d(match[i].x1, match[i].y1, -f);
            X2 = R * Eigen::Vector3d(match[i].x2, match[i].y2, -f);
            double N1 = (Bx * X2[2] - Bz * X2[0]) / (X1[0] * X2[2] - X2[0] * X1[2]);
            double N2 = (Bx * X1[2] - Bz * X1[0]) / (X1[0] * X2[2] - X2[0] * X1[2]);
            A(i, 0) = Bx;
            A(i, 1) = -X2[1] * Bx / X2[2];
            A(i, 2) = -X2[0] * X2[1] * N2 / X2[2];
            A(i, 3) = -(X2[2] + X2[1] * X2[1] / X2[2]) * N2;
            A(i, 4) = X2[0] * N2;
            L(i, 0) = N1 * X1[1] - N2 * X2[1] - By;
        }

        delta = (A.transpose() * A).inverse() * A.transpose() * L;


        u = u + delta[0];
        v = v + delta[1];
        phi = phi + delta[2];
        omega = omega + delta[3];
        kappa = kappa + delta[4];
        n++;
    } 
    while ((fabs(delta[0]) >= 3 * 10e-5 || fabs(delta[1]) >= 3 * 10e-5 || fabs(delta[2]) >= 3 * 10e-5 || fabs(delta[3]) >= 3 * 10e-5 || fabs(delta[4]) >= 3 * 10e-5) && (n < 60));
    cout << "迭代次数为:  " << n << endl;
    Eigen::VectorXd V(size);
    double vv = 0;
    double m0 = 0;
    V = A * delta - L;
    vv = V.dot(V);
    m0 = sqrt(vv / (size - 5));
    cout.precision(11);
    cout << "连续法相对定向求得相对定向元素为：" << endl
        << "u=" << u << endl
        << "v=" << v << endl
        << "phi=" << phi << endl
        << "omega=" << omega << endl
        << "kappa=" << kappa << endl;
    cout << "上下视差中误差为：" << endl
        << "m0=" << m0 << endl;
}