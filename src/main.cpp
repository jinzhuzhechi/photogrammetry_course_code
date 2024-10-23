#include <iostream>
#include <string>
#include <vector>
#include <fstream>


#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"

#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"

#include "photo.h"
#include "createEpiImg.h"
#include "relative_orientation.h"
#include "task.h"

#define PI 3.1415926

using namespace std;
using namespace cv;

int main(int, char**)
{
   int task = 10;
   // 后方交会
   if(task==1)
      return task1();
   // 前方交汇
   if(task==2)
      return task2();
   // 像片内定向 
   if(task==3)
      return task3();
   // 相对定向
   if(task==4)
      return task4();
   // 绝对定向 
   if(task==5)
      return task5();
   // 核线影像重构
   if(task==6)
      return task6();
   // 连续法相对定向
   if(task==7)
      return task7();
   // 空间后方交会，包含畸变参数
   if(task==8)
      return task8();
   // 直接线性变换
   if(task==9)
      return task9();
   // 相关系数匹配和最小二乘匹配
   if(task==10)
      return task10();
   return 0;
}

