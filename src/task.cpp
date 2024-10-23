#include "task.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"

#include "photo.h"
#include "createEpiImg.h"
#include "relative_orientation.h"
#include "data_reader.h"
#include "photo_pair.h"
#include "space_resection.h"
#include "dlt.h"
#include "matcher.h"
#include "detector.h"

using namespace std;
using namespace cv;

#define pi 3.1415926

//后方交会
int task1(void)
{
    //  初始化一张相片
    Photo aphoto = Photo(39800, 27500, 7600);

    string path_ctl_photo = "../data/1-后方交会/data1_ctl_photo.txt";
    string path_ctl_ground = "../data/1-后方交会/data1_ctl_ground.txt";
    int flag = datareader_1(path_ctl_photo, path_ctl_ground, aphoto);
    if (flag == -1)
        return -1;
    aphoto.resction();
    return 0;
}

//前方交汇
int task2()
{   
    Photo photo_left,photo_right;
    string path_ex = "../data/2-前方交会/extrinsics.txt"; 
    string path_in = "../data/2-前方交会/intrinsics.txt";  
    string path_l = "../data/2-前方交会/data2_321.txt"; 
    string path_r = "../data/2-前方交会/data2_320.txt"; 

    int flag = datareader_2(path_ex, path_in, path_l, path_r, photo_left, photo_right);
    if (flag == -1)
        return -1;

    //计算摄影基线
	double Bx = photo_right.X - photo_left.X;
	double By = photo_right.Y - photo_left.Y;
	double Bz = photo_right.Z - photo_left.Z;

    for (auto it = photo_left.ctl_photo.points.begin(); it != photo_left.ctl_photo.points.end(); ++it) 
    {
        int index = it->first;         // 获取键
        Mat_<double> pt1_3d_matrix = photo_left.getSRC(index);
        Mat_<double> pt2_3d_matrix = photo_right.getSRC(index);
        double N1_1 = (Bx * pt2_3d_matrix.at<double>(2, 0) - Bz * pt2_3d_matrix.at<double>(0, 0)) /
					(pt1_3d_matrix.at<double>(0, 0) * pt2_3d_matrix.at<double>(2, 0) - pt2_3d_matrix.at<double>(0, 0) * pt1_3d_matrix.at<double>(2, 0));
		double N2_1 = (Bx * pt1_3d_matrix.at<double>(2, 0) - Bz * pt1_3d_matrix.at<double>(0, 0)) /
					(pt1_3d_matrix.at<double>(0, 0) * pt2_3d_matrix.at<double>(2, 0) - pt2_3d_matrix.at<double>(0, 0) * pt1_3d_matrix.at<double>(2, 0));

        double N1 = (Bx * pt2_3d_matrix.at<double>(1, 0) - By * pt2_3d_matrix.at<double>(0, 0)) /
					(pt1_3d_matrix.at<double>(0, 0) * pt2_3d_matrix.at<double>(1, 0) - pt2_3d_matrix.at<double>(0, 0) * pt1_3d_matrix.at<double>(1, 0));
		double N2 = (Bx * pt1_3d_matrix.at<double>(2, 0) - By * pt1_3d_matrix.at<double>(0, 0)) /
					(pt1_3d_matrix.at<double>(0, 0) * pt2_3d_matrix.at<double>(1, 0) - pt2_3d_matrix.at<double>(0, 0) * pt1_3d_matrix.at<double>(1, 0));

		double Xt = photo_left.X + N1 * pt1_3d_matrix.at<double>(0, 0);
		double Yt = photo_left.Y + N1 * pt1_3d_matrix.at<double>(1, 0);
		double Zt = photo_left.Z + N1 * pt1_3d_matrix.at<double>(2, 0);
        cout<<index<<endl;
        cout<<"Xt: "<<Xt<<endl;
        cout<<"Yt: "<<Yt<<endl;
        cout<<"Zt: "<<Zt<<endl;
    }
    return 0;
}

// 像片内定向
int task3()
{
    struct FLAG {
		double x, y;
		double x_pixel, y_pixel;
		double x0, y0;
	}flag_info;

    // 输入基本参数
    double n_row = 11263; 
	double n_col = 11129;
	double pixel_size = 0.021;
	double x0 = 0.011, y0 = 0.002;

    ifstream file;
    file.open("../data/3-内定向/data3_320.txt");
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }

    string line;
    vector<vector<double>> data; // 存储读取的数据
    // 逐行读取文件
    while (std::getline(file, line)) 
    {
        istringstream iss(line);
        vector<double> row; // 存储当前行的数据
        // 从当前行的字符串流中提取数据
        double value;
        while (iss >> value) 
        {
            row.push_back(value);
        }
        // 将当前行的数据存入总体数据集
        data.push_back(row);
    }
    file.close();

    FLAG* flag = new FLAG[4];
    for (int i = 0; i < 4; i++)         //读入框标数据
	{
		flag[i].x = data[i][0]; 
        flag[i].y = data[i][1];
		flag[i].x_pixel = data[i][2]; 
        flag[i].y_pixel = data[i][3];
		flag[i].x0 = data[i][4];
        flag[i].y0 = data[i][5];
		//把单位从像素改成mm
		flag[i].x_pixel = (flag[i].x_pixel - n_col / 2) * pixel_size - x0;
		flag[i].y_pixel = (flag[i].y_pixel - n_row / 2) * pixel_size - y0;
	}

    //构建A矩阵
	Mat A = Mat::ones(1, 6, CV_64F);
	for (int i = 0; i < 4; i++)
	{
		double data[12] = { 1,flag[i].x_pixel,flag[i].y_pixel ,0,0,0,0,0,0,1,flag[i].x_pixel,flag[i].y_pixel };
		Mat temp = Mat(2, 6, CV_64F, data).clone();
		vconcat(A, temp, A);
	}
	A = A(cv::Range(1, 2 * 4 + 1), cv::Range::all());

	//构建l
	Mat l = Mat::ones(1, 1, CV_64F);
	for (int i = 0; i < 4; i++)
	{
		double data[2] = { flag[i].x,flag[i].y };
		Mat temp = Mat(2, 1, CV_64F, data).clone();
		vconcat(l, temp, l);
	}
	l = l(cv::Range(1, 2 * 4 + 1), cv::Range::all());

	//计算内定向变换参数
	Mat coff = (A.t() * A).inv() * (A.t() * l);

    //计算单位权中误差
	Mat v = A * coff - l;
	Mat temp = (v.t() * v);
	double result_sigma0 = sqrt(temp.at<double>(0, 0) / (2 * 4 - 6));
    cout<< result_sigma0<<endl;
    return 0;
}

// 独立模型法相对定向
int task4()
{
    Photo photo_left,photo_right;
    string path_ex = "../data/2-前方交会/extrinsics.txt"; 
    string path_in = "../data/2-前方交会/intrinsics.txt";  
    string path_l = "../data/2-前方交会/data2_321.txt"; 
    string path_r = "../data/2-前方交会/data2_320.txt"; 

    int flag = datareader_3(path_ex, path_in, path_l, path_r, photo_left, photo_right);
    if (flag == -1)
        return -1;
    photo_pair apair = photo_pair(photo_left,photo_right);
    apair.relative_orientation();
    return 0;
}

// 绝对定向
int task5()
{
    struct PT_PAIR {
		int ID;
		double Xp, Yp, Zp;
		double Xtp, Ytp, Ztp;
	};

    ifstream file;
    file.open("../data/5-绝对定向/data5.txt");
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
    string line;
    vector<vector<double>> data; // 存储读取的数据
    // 逐行读取文件
    while (std::getline(file, line)) 
    {
        istringstream iss(line);
        vector<double> row; // 存储当前行的数据
        // 从当前行的字符串流中提取数据
        double value;
        while (iss >> value) 
        {
            row.push_back(value);
        }
        // 将当前行的数据存入总体数据集
        data.push_back(row);
    }
    file.close();

    int n_pt = data.size();
    PT_PAIR *pt = new PT_PAIR[n_pt];
    for (int i = 0; i < data.size(); i++)
    {
        pt[i].ID = int(data[i][0]);
        pt[i].Xp = data[i][1];
        pt[i].Yp = data[i][2];
        pt[i].Zp = data[i][3];
        pt[i].Xtp = data[i][4];
        pt[i].Ytp = data[i][5];
        pt[i].Ztp = data[i][6];
    }
    	//给参数初值
	double lamda = 1;
	double phi = 0, omiga = 0, kappa = 0;
	double X0 = 0, Y0 = 0, Z0 = 0;

	double d1 = 1, d2 = 1, d3 = 1;  //判断角元素变化是否超过阈值
	int n_count = 0;   //迭代次数

	//计算重心化坐标
	double Xpg = 0, Ypg = 0, Zpg = 0;
	double Xtpg = 0, Ytpg = 0,Ztpg = 0;
	for (int i = 0; i < n_pt; i++)
	{
		Xpg += pt[i].Xp; Ypg += pt[i].Yp; Zpg += pt[i].Zp;
		Xtpg += pt[i].Xtp; Ytpg += pt[i].Ytp; Ztpg += pt[i].Ztp;
	}
	Xpg /= n_pt; Ypg /= n_pt; Zpg /= n_pt;
	Xtpg /= n_pt; Ytpg /= n_pt; Ztpg /= n_pt;

	for (int i = 0; i < n_pt; i++)
	{
		pt[i].Xp -= Xpg; pt[i].Yp -= Ypg; pt[i].Zp -= Zpg;
		pt[i].Xtp -= Xtpg; pt[i].Ytp -= Ytpg; pt[i].Ztp -= Ztpg;
	}

	while (d1 > 3e-5 || d2 > 3e-5 || d3 > 3e-5)
	{
		//构建误差方程
		Mat A = Mat::ones(1, 7, CV_64F);
		Mat l = Mat::ones(1, 1, CV_64F);
		double r[9];    //确定旋转矩阵

		for (int i = 0; i < n_pt; i++)
		{
		//更新旋转矩阵参数
			r[0] = cos(kappa) * cos(phi) - sin(phi) * sin(omiga) * sin(kappa);
			r[1] = -cos(phi) * sin(kappa) - sin(phi) * sin(omiga) * cos(kappa);
			r[2] = -sin(phi) * cos(omiga);
			r[3] = cos(omiga) * sin(kappa);
			r[4] = cos(kappa) * cos(omiga);
			r[5] = -sin(omiga);
			r[6] = sin(phi) * cos(kappa) + cos(phi) * sin(omiga) * sin(kappa);
			r[7] = -sin(phi) * sin(kappa) + cos(phi) * sin(omiga) * cos(kappa);
			r[8] = cos(phi) * cos(omiga);

			Mat rotate = Mat(3, 3, CV_64F, r).clone();
			double Xpg_array[3] = { pt[i].Xp ,pt[i].Yp ,pt[i].Zp };
			double Xtpg_array[3] = { pt[i].Xtp ,pt[i].Ytp ,pt[i].Ztp };
			Mat Xpg_matrix = Mat(3, 1, CV_64F, Xpg_array).clone();
			Mat Xtpg_matrix = Mat(3, 1, CV_64F, Xtpg_array).clone();

			Mat Xpg_matrix2 = rotate * Xpg_matrix;
			Mat Xpg_matrix3 = lamda * rotate * Xpg_matrix;

			double data_A[21] = { 1,0,0,Xpg_matrix2.at<double>(0,0) ,-lamda* Xpg_matrix2.at<double>(2,0) ,-Xpg_matrix3.at<double>(1,0)*sin(phi),-Xpg_matrix3.at<double>(1,0) * cos(phi) * cos(omiga)- Xpg_matrix3.at<double>(2,0) *sin(omiga),
			0,1,0,Xpg_matrix2.at<double>(1,0) ,0,Xpg_matrix3.at<double>(0,0)*sin(phi)- Xpg_matrix3.at<double>(2,0)*cos(phi),Xpg_matrix3.at<double>(0,0)* cos(phi)* cos(omiga)+ Xpg_matrix3.at<double>(2,0) * sin(phi) * cos(omiga),
			0,0,1,Xpg_matrix2.at<double>(2,0) ,lamda*Xpg_matrix2.at<double>(0,0) ,Xpg_matrix3.at<double>(1,0)* cos(phi),Xpg_matrix3.at<double>(0,0)*sin(omiga)- Xpg_matrix3.at<double>(1,0) * sin(phi) * cos(omiga) };
			Mat temp1 = Mat(3, 7, CV_64F, data_A).clone();
			vconcat(A, temp1, A);

			Mat temp2 = Xtpg_matrix - lamda * (rotate * Xpg_matrix);
			vconcat(l, temp2, l);
		}

		A = A(cv::Range(1, 3 * n_pt + 1), cv::Range::all());
		l = l(cv::Range(1, 3 * n_pt + 1), cv::Range::all());

		//平差并更新参数
		Mat delta = (A.t() * A).inv() * (A.t() * l);
		X0 = X0 + delta.at<double>(0, 0);
		Y0 = Y0 + delta.at<double>(1, 0);
		Z0 = Z0 + delta.at<double>(2, 0);
		lamda = lamda + delta.at<double>(3, 0);
		phi = phi + delta.at<double>(4, 0);
		omiga = omiga + delta.at<double>(5, 0);
		kappa = kappa + delta.at<double>(6, 0);

		d1 = abs(delta.at<double>(4, 0));
		d2 = abs(delta.at<double>(5, 0));
		d3 = abs(delta.at<double>(6, 0));
		n_count++;
	}
	int result_count = n_count;
	cout<<"迭代次数："<<result_count<<endl;
	cout<<"X0:"<<X0<<endl;
	cout<<"Y0:"<<Y0<<endl;
	cout<<"Z0:"<<Z0<<endl;
	cout<<"lamda:"<<lamda<<endl;
	cout<<"phi:"<<phi<<endl;
	cout<<"omiga:"<<omiga<<endl;
	cout<<"kappa:"<<kappa<<endl;
    return 0;
}




// 根据影像的外方位元素进行核线影像纠正
int task6(void)
{
    // string left_img_dir = "../data/aerial_simu/left.bmp";
   // string right_img_dir = "../data/aerial_simu/right.bmp";
   // IOP left_iop(512.00, 780.00, 50.0, 0.012, 1024, 1560);
   // IOP right_iop(512.00, 780.00, 50.0, 0.012, 1024, 1560);
   // EOP left_eop(2240.000000, 2002.000000, 8000.000000, 0.582000/180*PI, 1.210000/180*PI, 0.132000/180*PI);
   // EOP right_eop(2740.000000, 2002.000000, 8000.000000, -0.083000/180*PI, -0.311000/180*PI, 1.299000/180*PI);

   string left_img_dir = "../data/aerial_UAV/IMG_0182.jpg";
   string right_img_dir = "../data/aerial_UAV/IMG_0183.jpg";
   IOP left_iop(3744/2, 5616/2, 35.4904, 0.006400, 3744, 5616);
   IOP right_iop(3744 / 2, 5616 / 2, 35.4904, 0.006400, 3744, 5616);
   EOP left_eop(241608.3166, 3362629.0121, 551.2976, 0.057378, 0.088908, 0.046676);
   EOP right_eop(241693.5550, 3362627.8108, 551.2075, 0.054792, 0.084261, 0.003354);

   cv::Mat left_img = cv::imread(left_img_dir);
   cv::Mat right_img = cv::imread(right_img_dir);

   // Photo l_p(left_img, left_iop, left_eop);
   // Photo r_p(right_img, right_iop, right_eop);

   Photo epi_l;
   Photo epi_r;
   EOP epi_left_eop, epi_right_eop;

   // calculate_epiImg_EOP(l_p.eop, r_p.eop, epi_l.eop, epi_r.eop);
   calculate_epiImg_EOP(left_eop, right_eop, epi_left_eop, epi_right_eop);
   Mat left_M(3, 3, CV_64FC1);
   Mat right_M(3, 3, CV_64FC1);
   // calculate_M(left_eop, epi_l.eop, left_M);
   // calculate_M(right_eop, epi_r.eop, right_M);
   calculate_M(left_eop, epi_left_eop, left_M);
   calculate_M(right_eop, epi_right_eop, right_M);

   // 计算核线影像的尺寸
   int epi_width1, epi_height1, epi_width2, epi_height2, epi_width, epi_height;
   // calculate_epiImgSize(l_p.iop, left_M,  epi_width1, epi_height1);
   // calculate_epiImgSize(r_p.iop, right_M, epi_width2, epi_height2);
   calculate_epiImgSize(left_iop, left_M, epi_width1, epi_height1);
   calculate_epiImgSize(right_iop, right_M, epi_width2, epi_height2);

   epi_width = epi_width1 > epi_width2 ? epi_width1 : epi_width2;
   epi_height = epi_height1 > epi_height2 ? epi_height1 : epi_height2;

   // 制作核线影像
   //Mat left_epiImg(left_iop.height, left_iop.width, left_img.type());
   //Mat right_epiImg(right_iop.height, right_iop.width, right_img.type());
   Mat left_epiImg(epi_height, epi_width, left_img.type());
   Mat right_epiImg(epi_height, epi_width, right_img.type());
   creat_epiImg(left_img, left_epiImg, left_iop, left_M);
   creat_epiImg(right_img, right_epiImg, right_iop, right_M);

   int gap = 50;// 定义中间的空白区域宽度
   int height = std::max(left_epiImg.rows, right_epiImg.rows);  
   int width = left_epiImg.cols + right_epiImg.cols + gap;   
   cv::Mat combined_image = cv::Mat::zeros(height, width, left_epiImg.type());// 创建一个大小为 (height, width) 的新图像，并初始化为黑色
   left_epiImg.copyTo(combined_image(cv::Rect(0, 0, left_epiImg.cols, left_epiImg.rows)));
   right_epiImg.copyTo(combined_image(cv::Rect(left_epiImg.cols + gap, 0, right_epiImg.cols, right_epiImg.rows)));
   std::vector<int> lines_to_draw = { 20, 40, 60, 80 }; // 定义直线
   // 绘制直线
   for (int y : lines_to_draw) {
      int kk = y * combined_image.rows/ 100 ;
      cv::line(combined_image, cv::Point(0, kk), cv::Point(width, kk), cv::Scalar(0, 0, 255),3);
   }
   imwrite("combined_image.jpg", combined_image);
   crop_epiImg(combined_image);
   imwrite("combined_image_crop.jpg", combined_image);

   return 0;
}

// 连续法相对定向
int task7(void)
{
    string left_img_dir = "../data/aerial_UAV/IMG_0182.jpg";
    string right_img_dir = "../data/aerial_UAV/IMG_0183.jpg";
    cv::Mat left_img = cv::imread(left_img_dir);
    cv::Mat right_img = cv::imread(right_img_dir);
    get_relative_orientation(left_img, right_img, 35.4904);
    return 0;
}

// 空间后方交会，包含畸变参数
int task8(void)
{
    string gcp_path  = "../data/近景摄影测量/GCP.txt";
    string op_path_l = "../data/近景摄影测量/center_points_l.txt";
    string op_path_r = "../data/近景摄影测量/center_points_r.txt";
    string un_path_l = "../data/近景摄影测量/center_points_l_un.txt";
    string un_path_r = "../data/近景摄影测量/center_points_r_un.txt";

    map<int, GCP> gcps = readin_GCP(gcp_path); 
    map<int, GCP> gcps_l = readin_op(op_path_l, gcps);
    map<int, GCP> gcps_r = readin_op(op_path_r, gcps);
    map<int, UNP> un_p = readin_un(un_path_l, un_path_r);
    SpaceResection sr_l(4166, 0, 0, -2000, 0, -1000, 0, 0, 0);
    sr_l.gcps = gcps_l;
    sr_l.resection();

    SpaceResection sr_r(4166, 0, 0, -3600, 0, -1000, 0, 0, 0);
    sr_r.gcps = gcps_r;
    sr_r.resection();
    return 0;
}

// 直接线性变换
int task9(void)
{
    string gcp_path  = "../data/近景摄影测量/GCP.txt";
    string op_path_l = "../data/近景摄影测量/center_points_l.txt";
    string op_path_r = "../data/近景摄影测量/center_points_r.txt";
    string un_path_l = "../data/近景摄影测量/center_points_l_un.txt";
    string un_path_r = "../data/近景摄影测量/center_points_r_un.txt";

    map<int, GCP> gcps = readin_GCP(gcp_path); 
    map<int, GCP> gcps_l = readin_op(op_path_l, gcps);
    map<int, GCP> gcps_r = readin_op(op_path_r, gcps);
    map<int, UNP> un_p = readin_un(un_path_l, un_path_r);
    
    // 直接线性变换
    my_DLT dlt(gcps_l, gcps_r, un_p);
    dlt.dlt();

    return 0;
}

// 相关系数和最小二乘匹配
int task10(void)
{
    string result_corr_path = "result_corr.txt";
    string result_lsq_path = "result_lsq.txt";
    
    string path_l = "../data/match/2-town/l.jpg";
    string path_r = "../data/match/2-town/r.jpg";
    // string path_l = "../data/match/LOR49.bmp";
    // string path_r = "../data/match/LOR50.bmp";

// 读入图像
    Mat img_l = imread(path_l, IMREAD_COLOR);
    Mat img_r = imread(path_r, IMREAD_COLOR);

// 使用opencv中的各种函数直接得到两张影像的初始视差
    Mat homogeneous_matrix = Get_Homography_Matrix(img_l, img_r);
    int dx = static_cast<int>(std::round(homogeneous_matrix.at<double>(0, 2)));
    int dy = static_cast<int>(std::round(homogeneous_matrix.at<double>(1, 2)));
    cout<<homogeneous_matrix<<endl;

// 提取特征点，使用opencv的Maravec算法
    vector<cv::Point> corners_1, corners_2;
    vector<cv::Point> corners;

    // int corner_threshold = 700; // Moravec角点检测阈值
    int corner_threshold = 80; // Harris角点检测阈值

    // corners_1 = Moravec_Corner_Detector(img_l, 5, corner_threshold);
    // corners_2 = Moravec_Corner_Detector(img_r, 5, corner_threshold);

    corners_1 = Harris_Corner_Detector(img_l, corner_threshold);
    // corners_2 = Harris_Corner_Detector(img_r, corner_threshold);

    cout<<corners_1.size()<<endl;

    int corr_window_size = 25;    // 相关系数匹配窗口大小
    double corr_threshold = 0.95; // 相关系数匹配阈值

    int lsq_window_size = 5;     // 最小二乘匹配窗口大小
    double lsq_threshold = 0.95; //最小二乘匹配窗口大小

    // 显示角点
    Mat img_1_corner, img_2_corner;
    Draw_Corner(img_l, img_1_corner, corners_1);
    // Draw_Corner(img_r, img_2_corner, corners_2);

    cv::imwrite("角点1.jpg", img_1_corner);
    // cv::imwrite("角点2.jpg", img_2_corner);

    // 特征点匹配
    CorrelationMatcher corrMatcher;
    corrMatcher.setWindowSize(corr_window_size);
    corrMatcher.setThreshold(corr_threshold);

    vector<Match> corrMatches;
    corrMatcher.Correlation_Match_initiated(img_l, img_r, corners_1, dx, dy, corrMatches);
    // corrMatcher.Correlation_Match(img_l, img_r, corners_1, corners_2, corrMatches);
    // 显示匹配结果
    cv::Mat img_match;
    corrMatcher.drawMatches(img_l, img_r, img_match, corrMatches);
    cv::imshow("相关系数匹配结果", img_match);
    cv::waitKey(0);

    std::ofstream ofs;
    ofs.open(result_corr_path);
    for (auto match : corrMatches)
    {
        ofs << "srcX: " << match.srcPt.x
            << "\tsrcY: " << match.srcPt.y
            << "\tdstX: " << match.dstPt.x
            << "\tdstY: " << match.dstPt.y
            << "\tidx: " << match.dist << std::endl;
    }
    ofs.close();

    // 在相关系数匹配的基础上进行最小二乘匹配
    LsqMatcher lsqMatcher;
    lsqMatcher.setWindowSize(lsq_window_size);
    lsqMatcher.setThreshold(lsq_threshold);
    std::vector<Match> lsqMatches;

    for (auto match : corrMatches)
    {
        if (lsqMatcher.SubPixel_Match(img_l, img_r, match))
            lsqMatches.push_back(match);
    }

    img_match = cv::Mat();
    lsqMatcher.drawMatches(img_l, img_r, img_match, lsqMatches);
    cv::imshow("最小二乘匹配结果", img_match);
    cv::waitKey(0);
    ofs.open(result_lsq_path);
    for (auto match : lsqMatches)
    {
        ofs << "srcX: " << match.srcPt.x
            << "\tsrcY: " << match.srcPt.y
            << "\tdstX: " << match.dstPt.x
            << "\tdstY: " << match.dstPt.y
            << "\tidx: " << match.dist << std::endl;
    }
    ofs.close();
    return 0;
}

