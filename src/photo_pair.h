#ifndef PHOTO_PAIR
#define PHOTO_PAIR

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"

#include "photo.h"
using namespace cv;
using namespace std;


class photo_pair
{
public:
    Photo photo_left;
    Photo photo_right;

    double u = 0, v = 0;
	double Phi = 0, Omega = 0, Kappa =0;
    double Bx,By,Bz;

    // 定义像对间的旋转矩阵
    Mat tran;
    vector<int> index;

    vector<int> get_ctl_index()
    {
        vector<int> keys;
        for (const auto& pair : photo_left.ctl_photo.points) 
        {
            if (photo_right.ctl_photo.points.count(pair.first) > 0) 
            {
                int key = pair.first;
                keys.push_back(key);
            }
        }
        return keys;
    }

    photo_pair() = default;
    photo_pair(Photo photo1, Photo photo2)
    {
        this->photo_left = photo1;
        this->photo_right = photo2;
        index = get_ctl_index();
        this->tran = transform();
    }
    void relative_orientation();

private:
    Mat transform();
};



#endif