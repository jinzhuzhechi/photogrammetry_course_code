#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"

#include "photo.h"

using namespace cv;
using namespace std;

// 对相片对象进行初始化
Photo::Photo(double x, double y, double z, double phi, double omega, double kappa, double x_0_in, double y_0_in, double f_in, int m_in)
{
    this->X = x;
    this->Y = y;
    this->Z = z;
    this->Phi = phi;
    this->Omega = omega;
    this->Kappa = kappa;
    this->x_0 = x_0_in;
    this->y_0 = y_0_in;
    this->f = f_in;
    this->m = m_in;
    this->tran = transform();
}

// 获取旋转矩阵
Mat Photo::transform()
{
    Mat tran = Mat::zeros(3,3,CV_64FC1);
    tran.at<double>(0,0) = cos(this->Phi)*cos(this->Kappa)-sin(this->Phi)*sin(this->Omega)*sin(this->Kappa);
    tran.at<double>(0,1) = -cos(this->Phi)*sin(this->Kappa)-sin(this->Phi)*sin(this->Omega)*cos(this->Kappa);
    tran.at<double>(0,2) = -sin(this->Phi)*cos(this->Omega);

    tran.at<double>(1,0) = cos(this->Omega)*sin(this->Kappa);
    tran.at<double>(1,1) = cos(this->Omega)*cos(this->Kappa);
    tran.at<double>(1,2) = -sin(this->Omega);

    tran.at<double>(2,0) = sin(this->Phi)*cos(this->Kappa)+cos(this->Phi)*sin(this->Omega)*sin(this->Kappa);
    tran.at<double>(2,1) = -sin(this->Phi)*sin(this->Kappa)+cos(this->Phi)*sin(this->Omega)*cos(this->Kappa);
    tran.at<double>(2,2) = cos(this->Phi)*cos(this->Omega);
    return tran;
};


// 进行后方交会
void Photo::resction()
{
    Mat_<double>delta_Extrinsics(6, 1);
    do
    {
        Mat_<double>p_xy(ctl_ground.getNumPoints(), 2);
        Mat_<double>XYZ_bar(3,ctl_ground.getNumPoints());
        for(int i=0;i<ctl_ground.getNumPoints();i++)
        {
            Mat_<double>p_delta_XYZ(3, 1);
            Mat_<double> p_center(3, 1);
            p_center(0, 0) = X;
            p_center(1, 0) = Y;
            p_center(2, 0) = Z;
            p_delta_XYZ = ctl_ground.getPoint(i) - p_center;
            Mat_<double> inverse_tran;
            invert(this->tran, inverse_tran);
            XYZ_bar.col(i) = inverse_tran * p_delta_XYZ;
            p_xy(i,0) = - this->f * XYZ_bar.at<double>(0,i) / XYZ_bar.at<double>(2,i);
            p_xy(i,1) = - this->f * XYZ_bar.at<double>(1,i) / XYZ_bar.at<double>(2,i);
        }
        Mat_<double>A(ctl_ground.getNumPoints()*2, 6);
        for(int i=0;i<ctl_ground.getNumPoints();i++)
        {
            A(2*i,0) = (this->tran.at<double>(0,0)*this->f + this->tran.at<double>(0,2)*p_xy.at<double>(i,0)) / XYZ_bar.at<double>(2,i);
            A(2*i,1) = (this->tran.at<double>(1,0)*this->f + this->tran.at<double>(1,2)*p_xy.at<double>(i,0)) / XYZ_bar.at<double>(2,i);
            A(2*i,2) = (this->tran.at<double>(2,0)*this->f + this->tran.at<double>(2,2)*p_xy.at<double>(i,0)) / XYZ_bar.at<double>(2,i);
            A(2*i,3) = p_xy.at<double>(i,1)*sin(this->Omega)-((p_xy.at<double>(i,0)*cos(this->Kappa) - p_xy.at<double>(i,1)*sin(this->Kappa))*p_xy.at<double>(i,0)/this->f + this->f*cos(this->Kappa))*cos(this->Omega);
            A(2*i,4) = -this->f*sin(this->Kappa) - (p_xy.at<double>(i,0)*sin(this->Kappa) + p_xy.at<double>(i,1)*cos(this->Kappa))*p_xy.at<double>(i,0)/this->f;
            A(2*i,5) = p_xy.at<double>(i,1);

            A(2*i+1,0) = (this->tran.at<double>(0,1)*this->f + this->tran.at<double>(0,2)*p_xy.at<double>(i,1)) / XYZ_bar.at<double>(2,i);
            A(2*i+1,1) = (this->tran.at<double>(1,1)*this->f + this->tran.at<double>(1,2)*p_xy.at<double>(i,1)) / XYZ_bar.at<double>(2,i);
            A(2*i+1,2) = (this->tran.at<double>(2,1)*this->f + this->tran.at<double>(2,2)*p_xy.at<double>(i,1)) / XYZ_bar.at<double>(2,i);
            // A(2*i+1,3) = -p_xy.at<double>(i,0)*sin(this->Omega)-((p_xy.at<double>(i,0)*cos(this->Kappa) - p_xy.at<double>(i,1)*sin(this->Kappa))*p_xy.at<double>(i,1)/this->f - this->f*cos(this->Kappa))*cos(this->Omega);
            // A(2*i+1,4) = -this->f*cos(this->Kappa) - (p_xy.at<double>(i,0)*sin(this->Kappa) + p_xy.at<double>(i,1)*cos(this->Kappa))*p_xy.at<double>(i,1)/this->f;
            A(2*i+1,3) = -(p_xy.at<double>(i,0) * p_xy.at<double>(i,1))*tran.at<double>(1,1) / f + (f + pow(p_xy.at<double>(i,1), 2) / f)*tran.at<double>(1,0) - p_xy.at<double>(i,0) * tran.at<double>(1,2);
            A(2*i+1,4) = -pow(p_xy.at<double>(i,1), 2)*cos(Kappa) / f - (p_xy.at<double>(i,0) * p_xy.at<double>(i,1)) / f*cos(Kappa) - f*sin(Kappa);
            A(2*i+1,5) = -p_xy.at<double>(i,0);
        }
        Mat_<double>L(ctl_ground.getNumPoints()*2, 1);
        for(int i=0;i<ctl_ground.getNumPoints()*2;i++)
        {
            if (i % 2 == 0)	//x
                L(i,0) = ctl_photo.getPoint(i / 2).at<double>(0, 0)- p_xy.at<double>(i / 2,0);
            else
                L(i,0) = ctl_photo.getPoint(i / 2).at<double>(1, 0)- p_xy.at<double>(i / 2,1);
        }

        Mat_<double>ATA = A.t() * A;
        Mat_<double>ATL = A.t() * L;
        Mat_<double>inv_ATA = ATA.inv();
        delta_Extrinsics = inv_ATA * ATL;
        renewExtrinsics(delta_Extrinsics);
    }
    while (delta_Extrinsics.at<double>(3,0)>3e-5 || delta_Extrinsics.at<double>(4,0)>3e-5 || delta_Extrinsics.at<double>(5,0)>3e-5);

    showExtrinsics();
}

void Photo::renewExtrinsics(Mat_<double> delta_extrinsics)
{
    this->X = this->X+delta_extrinsics.at<double>(0,0);
    this->Y = this->Y+delta_extrinsics.at<double>(1,0);
    this->Z = this->Z+delta_extrinsics.at<double>(2,0);
    this->Phi = this->Phi+delta_extrinsics.at<double>(3,0);
    this->Omega = this->Omega+delta_extrinsics.at<double>(4,0);
    this->Kappa = this->Kappa+delta_extrinsics.at<double>(5,0);
    this->tran = transform();
}

void Photo::showExtrinsics()
{
    cout<<"X: "<<this->X<<endl;
    cout<<"Y: "<<this->Y<<endl;
    cout<<"Z: "<<this->Z<<endl;
    cout<<"tran: "<<this->tran<<endl;
}

// 获取像空间辅助坐标Spatial reference coordinates
Mat_<double> Photo::getSRC(int index)
{
    Mat_<double> xyz(3, 1);
    xyz(0, 0) = ctl_photo.getPoint(index).at<double>(0, 0) - this->x_0;
    xyz(1, 0) = ctl_photo.getPoint(index).at<double>(1, 0) - this->y_0;
    xyz(2, 0) = -this->f;
    return tran * xyz;
}



