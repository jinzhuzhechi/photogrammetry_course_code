#include "createEpiImg.h"

//计算核线影像的外方位元素
void calculate_epiImg_EOP(EOP left_eop, EOP right_eop, EOP& left_epiImg_eop, EOP& right_epiImg_eop) 
{
    double BX = right_eop.Xs - left_eop.Xs;
    double BY = right_eop.Ys - left_eop.Ys;
    double BZ = right_eop.Zs - left_eop.Zs;
    
    double B = sqrt(BX * BX + BY * BY + BZ * BZ);
    double T = atan2(BY, BX);
    double V = asin(BZ / B);
    
    double sint = sin(T), cost = cos(T);
    double sinv = sin(V), cosv = cos(V);

    right_epiImg_eop.a1 = left_epiImg_eop.a1 = cost * cosv;
    right_epiImg_eop.a2 = left_epiImg_eop.a2 = sint * cosv;
    right_epiImg_eop.a3 = left_epiImg_eop.a3 = sinv;
    right_epiImg_eop.b1 = left_epiImg_eop.b1 = -sint;
    right_epiImg_eop.b2 = left_epiImg_eop.b2 = cost;
    right_epiImg_eop.b3 = left_epiImg_eop.b3 = 0;
    right_epiImg_eop.c1 = left_epiImg_eop.c1 = -cost * sinv;
    right_epiImg_eop.c2 = left_epiImg_eop.c2 = -sint * sinv;
    right_epiImg_eop.c3 = left_epiImg_eop.c3 = cosv;
}

//计算旋转矩阵M
void calculate_M(EOP eop, EOP epiImg_eop, Mat& M) 
{
    Mat AA = (Mat_<double>(3, 3) << epiImg_eop.a1, epiImg_eop.a2, epiImg_eop.a3, epiImg_eop.b1, epiImg_eop.b2, epiImg_eop.b3, epiImg_eop.c1, epiImg_eop.c2, epiImg_eop.c3);
    Mat BB = (Mat_<double>(3, 3) << eop.a1, eop.a2, eop.a3, eop.b1, eop.b2, eop.b3, eop.c1, eop.c2, eop.c3);
    Mat CC;
    transpose(AA, CC);
    Mat DD;
    transpose(BB, DD);
    M = DD * CC;
    // 输出结果
    // cout << "M: " << M << endl;
    // cout << AA.transpose() << endl;
}


// 计算核线影像的尺寸
void calculate_epiImgSize(IOP iop, Mat left_M, int& epi_width, int& epi_height) {
    // 获取原始影像内方位元素
    int nWidth = iop.width;
    int nHeight = iop.height;
    double f = iop.f;
    double fn = iop.f;
    double pixelSize = iop.pixelSize;
    // 计算原始图像图像中心点坐标
    double x0 = nWidth / 2;
    double y0 = nHeight / 2;
    // 计算原始影像四个角点坐标
    double ori_UL_x = -x0 * pixelSize, ori_UL_y = y0 * pixelSize;
    double ori_UR_x = x0 * pixelSize, ori_UR_y = y0 * pixelSize;
    double ori_DL_x = -x0 * pixelSize, ori_DL_y = -y0 * pixelSize;
    double ori_DR_x = x0 * pixelSize, ori_DR_y = -y0 * pixelSize;
    // 计算旋转后影像四个角点坐标
    double UL_x = -fn * (left_M.at<double>(0, 0) * ori_UL_x + left_M.at<double>(0, 1) * ori_UL_y - left_M.at<double>(0, 2) * f) / (left_M.at<double>(2, 0) * ori_UL_x + left_M.at<double>(2, 1) * ori_UL_y - left_M.at<double>(2, 2) * f);
    double UL_y = -fn * (left_M.at<double>(1, 0) * ori_UL_x + left_M.at<double>(1, 1) * ori_UL_y - left_M.at<double>(1, 2) * f) / (left_M.at<double>(2, 0) * ori_UL_x + left_M.at<double>(2, 1) * ori_UL_y - left_M.at<double>(2, 2) * f);
    double UR_x = -fn * (left_M.at<double>(0, 0) * ori_UR_x + left_M.at<double>(0, 1) * ori_UR_y - left_M.at<double>(0, 2) * f) / (left_M.at<double>(2, 0) * ori_UR_x + left_M.at<double>(2, 1) * ori_UR_y - left_M.at<double>(2, 2) * f);
    double UR_y = -fn * (left_M.at<double>(1, 0) * ori_UR_x + left_M.at<double>(1, 1) * ori_UR_y - left_M.at<double>(1, 2) * f) / (left_M.at<double>(2, 0) * ori_UR_x + left_M.at<double>(2, 1) * ori_UR_y - left_M.at<double>(2, 2) * f);
    double DL_x = -fn * (left_M.at<double>(0, 0) * ori_DL_x + left_M.at<double>(0, 1) * ori_DL_y - left_M.at<double>(0, 2) * f) / (left_M.at<double>(2, 0) * ori_DL_x + left_M.at<double>(2, 1) * ori_DL_y - left_M.at<double>(2, 2) * f);
    double DL_y = -fn * (left_M.at<double>(1, 0) * ori_DL_x + left_M.at<double>(1, 1) * ori_DL_y - left_M.at<double>(1, 2) * f) / (left_M.at<double>(2, 0) * ori_DL_x + left_M.at<double>(2, 1) * ori_DL_y - left_M.at<double>(2, 2) * f);
    double DR_x = -fn * (left_M.at<double>(0, 0) * ori_DR_x + left_M.at<double>(0, 1) * ori_DR_y - left_M.at<double>(0, 2) * f) / (left_M.at<double>(2, 0) * ori_DR_x + left_M.at<double>(2, 1) * ori_DR_y - left_M.at<double>(2, 2) * f);
    double DR_y = -fn * (left_M.at<double>(1, 0) * ori_DR_x + left_M.at<double>(1, 1) * ori_DR_y - left_M.at<double>(1, 2) * f) / (left_M.at<double>(2, 0) * ori_DR_x + left_M.at<double>(2, 1) * ori_DR_y - left_M.at<double>(2, 2) * f);

    // 确定核线影像范围
    double max_x = max(max(UL_x, UR_x), max(DL_x, DR_x));
    double max_y = max(max(UL_y, UR_y), max(DL_y, DR_y));
    double min_x = min(min(UL_x, UR_x), min(DL_x, DR_x));
    double min_y = min(min(UL_y, UR_y), min(DL_y, DR_y));

    epi_width = max_x > -min_x ? (2 * max_x / pixelSize) : (-2 * min_x / pixelSize);
    epi_height = max_y > -min_y ? (2 * max_y / pixelSize) : (-2 * min_y / pixelSize);
}



// 制作核线影像
void creat_epiImg(Mat img, Mat& epiImg, IOP iop, Mat M) 
{

    double u, v, x, y;
    int c, r;
    double x1, x2, y1, y2;
    double f = iop.f;
    double pixelSize = iop.pixelSize;
    double xo = iop.cx;
    double yo = iop.cy;

    for (int j = 0; j < epiImg.rows; j++) {//epi_height
        for (int i = 0; i < epiImg.cols; i++) {//epi_width
            u = i * pixelSize - epiImg.cols * pixelSize / 2;//核线影像坐标
            v = j * pixelSize - epiImg.rows * pixelSize / 2;
            x = -f * (M.at<double>(0, 0) * u + M.at<double>(1, 0) * v - M.at<double>(2, 0) * f ) / (M.at<double>(0, 2) * u + M.at<double>(1, 2) * v - M.at<double>(2, 2) * f );//倾斜影像坐标
            y = -f * (M.at<double>(0, 1) * u + M.at<double>(1, 1) * v - M.at<double>(2, 1) * f ) / (M.at<double>(0, 2) * u + M.at<double>(1, 2) * v - M.at<double>(2, 2) * f );
            c = x / pixelSize + xo;
            r = y / pixelSize + yo;

            if ((c < 0.0) || (r < 0.0) || (r >= iop.height - 1) || (c >= iop.width - 1)) {
                epiImg.at<cv::Vec3b>(j, i)[0] = 0;
                epiImg.at<cv::Vec3b>(j, i)[1] = 0;
                epiImg.at<cv::Vec3b>(j, i)[2] = 0;
            }
            else {//进行双线性内插：f(x, y) = f(0, 0)(1 - x)(1 - y) + f(1, 0)x(1 - y) + f(0, 1)(1 - x)y + f(1, 1)xy
                x1 = floor(c), x2 = ceil(c), y1 = floor(r), y2 = ceil(r);
                if (x1 == x2 && y1 != y2) { // x1 和 x2 相同，y1 和 y2 不同
                    for (int k = 0; k < 3; ++k) { // 对每个通道进行插值
                        epiImg.at<cv::Vec3b>(j, i)[k] = img.at<cv::Vec3b>(y1, x1)[k] * (y2 - r) + img.at<cv::Vec3b>(y2, x1)[k] * (r - y1);
                    }
                    
                }
                else if (x1 != x2 && y1 == y2) { // x1 和 x2 相同，y1 和 y2 相同
                    for (int k = 0; k < 3; ++k) { // 对每个通道进行赋值
                        epiImg.at<cv::Vec3b>(j, i)[k] = img.at<cv::Vec3b>(y1, x1)[k] * (x2 - c) + img.at<cv::Vec3b>(y1, x2)[k] * (c - x1);
                    }
                }
                else if (x1 == x2 && y1 == y2) { // x1 和 x2 相同，y1 和 y2 相同
                    for (int k = 0; k < 3; ++k) { // 对每个通道进行赋值
                        epiImg.at<cv::Vec3b>(j, i)[k] = img.at<cv::Vec3b>(y1, x1)[k];
                    }
                }
                else if (x1 != x2 && y1 != y2) { // x1 和 x2 不同，y1 和 y2 不同
                    for (int k = 0; k < 3; ++k) { // 对每个通道进行插值
                        epiImg.at<cv::Vec3b>(j, i)[k] = img.at<cv::Vec3b>(y1, x1)[k] * (x2 - c) * (y2 - r) + img.at<cv::Vec3b>(y1, x2)[k] * (c - x1) * (y2 - r) +
                            img.at<cv::Vec3b>(y2, x1)[k] * (x2 - c) * (r - y1) + img.at<cv::Vec3b>(y2, x2)[k] * (c - x1) * (r - y1);
                    }
                }
            }
        }
    }
}


void crop_epiImg(Mat& image) {
    Mat gray;
    cvtColor(image, gray, COLOR_BGR2GRAY);
    Mat thresh;
    threshold(gray, thresh, 0, 255,THRESH_OTSU);    
    std::vector<Point> nonzero_points;
    findNonZero(thresh, nonzero_points);// 找到图像中非零像素点的坐标   
    Rect bounding_rect = boundingRect(nonzero_points);// 找到最小的包围矩形
    image = image(bounding_rect);  
}