#ifndef TASK
#define TASK

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

// 后方交会
int task1(void);
// 前方交汇
int task2(void);
// 像片内定向
int task3(void);
// 相对定向
int task4(void);
// 绝对定向
int task5(void);
// 核线影像重构
int task6(void);
// 连续法相对定向
int task7(void);
// 空间后方交会，包含畸变参数
int task8(void);
// 直接线性变换
int task9(void);
// 相关系数匹配和最小二乘匹配
int task10(void);




#endif