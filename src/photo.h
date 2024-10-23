#ifndef PHOTO
#define PHOTO

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"

#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"

using namespace cv;
using namespace std;


// 控制点对象
class ControlPoints 
{
public:
    ControlPoints() = default;
    void addPoint(int index, Mat_<double> newPoint) 
    {
        points[index] = newPoint;
    }
    
    void removePoint(int index) 
    {
        points.erase(index);
    }
    
    int getNumPoints() const 
    {
        return points.size();
    }
    
    Mat_<double> getPoint(int index) const 
    {
        auto it = points.find(index);
        if (it != points.end()) 
        {
            return it->second;
        } 
        else 
        {
            // 返回一个默认的 my_point 对象，可以根据实际需求进行修改
            return Mat_<double>();
        }
    }
    std::unordered_map<int, Mat_<double> > points;
};


// 内定向元素
struct IOP {
    double cx;//像平面坐标原点对应的扫描坐标值
    double cy;
    double plane2scan[4];//像平面坐标到扫描坐标转换系数
    double f;//主距
    double pixelSize;//像素大小
    int width;
    int height;

    IOP(){}
    IOP(double cx_, double cy_, double f_, double pixelSize_, int width_, int height_) :cx(cx_), cy(cy_), f(f_), pixelSize(pixelSize_), width(width_), height(height_) {
        plane2scan[0] = plane2scan[3] = 1;
        plane2scan[1] = plane2scan[2] = 0;
    }
    IOP(double cx_, double cy_, double f_) :cx(cx_), cy(cy_), f(f_) {
        plane2scan[0] = plane2scan[3] = 1;
        plane2scan[1] = plane2scan[2] = 0;
    }
};

// 外方位元素
struct EOP
{
    double Xs, Ys, Zs;  // 图像的空间位置
    double Phi, Omega, Kappa;  // 旋转角度
    double a1, a2, a3;
    double b1, b2, b3;
    double c1, c2, c3;
    Mat tran;
    EOP() {}
    EOP(double Xs_, double Ys_, double Zs_, double Phi_, double Omega_, double Kappa_) :Xs(Xs_), Ys(Ys_), Zs(Zs_), Phi(Phi_), Omega(Omega_), Kappa(Kappa_) 
    {
        double cosp = cos(Phi);    double sinp = sin(Phi);
        double cosw = cos(Omega);  double sinw = sin(Omega);
        double cosk = cos(Kappa);  double sink = sin(Kappa);
        a1 = cosp * cosk - sinp * sinw * sink;
        a2 = -cosp * sink - sinp * sinw * cosk;
        a3 = -sinp * cosw;
        b1 = cosw * sink;
        b2 = cosw * cosk;
        b3 = -sinw;
        c1 = sinp * cosk + cosp * sinw * sink;
        c2 = -sinp * sink + cosp * sinw * cosk;
        c3 = cosp * cosw;
        tran = (Mat_<double>(3, 3) << a1, a2, a3, b1, b2, b3, c1, c2, c3);
    }
};

class Photo
{
public:
// 定义外部参数
    double X, Y, Z, Phi,Omega,Kappa;
// 定义内部参数
    double x_0, y_0, f;
// 比例尺
    int m;
    // 航带号
    int band = 0;
    // 像片号
    int index = 0;
    // 定义从像平面到像空间辅助坐标系的辅助矩阵
    Mat tran;
    // 定义像控制点和其对应地面点坐标
    ControlPoints ctl_photo;
    ControlPoints ctl_ground;
    // 获取像空间辅助坐标Spatial reference coordinates
    Mat_<double> getSRC(int index);

    Photo() = default;
    Photo(double x, double y, double z, double phi = 0, double omega = 0, double kappa = 0, double x_0_in = 0, double y_0_in = 0 , double f_in = 0.15324, int m_in = 40000);
    ~Photo(){};
    void resction();
    void renewExtrinsics(Mat_<double> extrinsics);
    void showExtrinsics();

private:
    Mat transform();
};

#endif
