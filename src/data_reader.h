#ifndef DATA_READER
#define DATA_READER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"

using namespace std;
using namespace cv;

#include "photo.h"

struct GCP
{
    double x;
    double y;
    double X;
    double Y;
    double Z;
};

struct UNP
{
    double x_l;
    double y_l;
    double x_r;
    double y_r;
    double X;
    double Y;
    double Z;
};



// 前方交会的数据读取
int datareader_1(string path_ctl_photo, string path_ctl_ground, Photo& aphoto);
int datareader_2(string path_ex, string path_in, string path_l, string path_r, Photo& photo_left, Photo& photo_right);
int datareader_3(string path_ex, string path_in, string path_l, string path_r, Photo& photo_left, Photo& photo_right);

map<int, GCP> readin_GCP(string file_path);
// 读入observed_point
map<int, GCP> readin_op(string file_path, map<int, GCP> gcps);

map<int, UNP> readin_un(string file_path_l, string file_path_r);


#endif