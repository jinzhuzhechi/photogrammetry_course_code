
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/core.hpp"

#include "photo.h"
#include "photo_pair.h"

using namespace cv;
using namespace std;

// 获取旋转矩阵
Mat photo_pair::transform()
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

void photo_pair:: relative_orientation()
{
    Mat_<double> A(index.size(),5);
    Mat_<double> L(index.size(),1);
    Mat_<double> delta_pram(5,1);
    do
    {
        // 基线暂定为1
        Bx = 1;
        By = By * u; 
        Bz = By * v;
        tran = transform();
        for(int i = 0;i<index.size();i++)
        {
            Mat_<double> point_left_xy = photo_left.ctl_photo.getPoint(index[i]);
            Mat_<double> point_right_xy = photo_right.ctl_photo.getPoint(index[i]);

            Mat_<double> point_left_xyf(3,1);
            point_left_xyf.at<double>(0,0) = point_left_xy.at<double>(0,0)-photo_left.x_0;
            point_left_xyf.at<double>(1,0) = point_left_xy.at<double>(1,0)-photo_left.y_0;
            point_left_xyf.at<double>(2,0) = -photo_left.f;
            Mat_<double> point_right_xyf(3,1);
            point_right_xyf.at<double>(0,0) = point_right_xy.at<double>(0,0)-photo_right.x_0;
            point_right_xyf.at<double>(1,0) = point_right_xy.at<double>(1,0)-photo_right.y_0;
            point_right_xyf.at<double>(2,0) = -photo_right.f;

            Mat_<double> p_l_XYZ = point_left_xyf;
            Mat_<double> p_r_XYZ = tran * point_right_xyf;

            //计算点投影系数
            double N1 = (Bx * p_r_XYZ.at<double>(2, 0) - Bz * p_r_XYZ.at<double>(0, 0)) /
                (p_l_XYZ.at<double>(0, 0) * p_r_XYZ.at<double>(2, 0) - p_r_XYZ.at<double>(0, 0) * p_l_XYZ.at<double>(2, 0));
            double N2 = (Bx * p_l_XYZ.at<double>(2, 0) - Bz * p_l_XYZ.at<double>(0, 0)) /
                (p_l_XYZ.at<double>(0, 0) * p_r_XYZ.at<double>(2, 0) - p_r_XYZ.at<double>(0, 0) * p_l_XYZ.at<double>(2, 0));

            // double N1 = (Bx * p_r_XYZ.at<double>(1, 0) - By * p_r_XYZ.at<double>(0, 0)) /
            //     (p_l_XYZ.at<double>(0, 0) * p_r_XYZ.at<double>(1, 0) - p_r_XYZ.at<double>(0, 0) * p_l_XYZ.at<double>(1, 0));
            // double N2 = (Bx * p_l_XYZ.at<double>(1, 0) - By * p_l_XYZ.at<double>(0, 0)) /
            //     (p_l_XYZ.at<double>(0, 0) * p_r_XYZ.at<double>(1, 0) - p_r_XYZ.at<double>(0, 0) * p_l_XYZ.at<double>(1, 0));

            A.at<double>(i,0) = Bx;
            A.at<double>(i,1) = -Bx * p_r_XYZ.at<double>(1,0) / p_r_XYZ.at<double>(2,0);
            A.at<double>(i,2) = -N2 * p_r_XYZ.at<double>(0,0) * p_r_XYZ.at<double>(1,0) / p_r_XYZ.at<double>(2,0);
            A.at<double>(i,3) = -N2 * (p_r_XYZ.at<double>(2,0) + p_r_XYZ.at<double>(1,0) * p_r_XYZ.at<double>(1,0) / p_r_XYZ.at<double>(2,0));
            A.at<double>(i,4) = p_r_XYZ.at<double>(0,0) * N2;
            L.at<double>(i,0) = N1 * p_l_XYZ.at<double>(1, 0) - N2 * p_r_XYZ.at<double>(1, 0) - By;

            cout<<"N1: "<<N1<<endl;
            cout<<"N2: "<<N2<<endl;
        }

        delta_pram = (A.t() * A).inv() * (A.t() * L);
        u = u + delta_pram.at<double>(0,0);
        v = v + delta_pram.at<double>(1,0);
        Phi = Phi + delta_pram.at<double>(2,0);
        Omega = Omega + delta_pram.at<double>(3,0);
        Kappa = Kappa + delta_pram.at<double>(4,0);
    } 
    while(abs(delta_pram.at<double>(2,0)) > 3e-5 || abs(delta_pram.at<double>(3,0)) > 3e-5 || abs(delta_pram.at<double>(4,0)) > 3e-5);
    cout<<"u: "<<u<<endl;
    cout<<"v: "<<v<<endl;
    cout<<"Phi: "<<Phi<<endl;
    cout<<"Omega: "<<Omega<<endl;
    cout<<"Kappa: "<<Kappa<<endl;
}