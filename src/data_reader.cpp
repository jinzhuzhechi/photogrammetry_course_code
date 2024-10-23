#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"

using namespace std;
using namespace cv;

#include "data_reader.h"
#include "photo.h"

#define pi 3.1415926

int datareader_1(string path_ctl_photo, string path_ctl_ground, Photo& aphoto)
{
    ifstream file;
    file.open(path_ctl_photo);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
    string line;
    getline(file, line);
    int num_point = stoi(line);

    for (int i = 0; i < num_point; i++)
    {
        getline(file, line);
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (getline(ss, cell, ','))
        {
            row.push_back(cell);
        }
        Mat point(2, 1, CV_64F);
        point.at<double>(0, 0) = stod(row[0]);
        point.at<double>(1, 0) = stod(row[1]);

        aphoto.ctl_photo.addPoint(i, point);
    }
    file.close();

    file.open(path_ctl_ground);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
    getline(file, line);
    num_point = stoi(line);
    for (int i = 0; i < num_point; i++)
    {
        getline(file, line);
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> row;
        while (getline(ss, cell, ','))
        {
            row.push_back(cell);
        }
        Mat point(3, 1, CV_64F);
        point.at<double>(0, 0) = stod(row[0]);
        point.at<double>(1, 0) = stod(row[1]);
        point.at<double>(2, 0) = stod(row[2]);

        aphoto.ctl_ground.addPoint(i, point);
    }
    return 0;
}


int datareader_2(string path_ex, string path_in, string path_l, string path_r, Photo& photo_left, Photo& photo_right)
{
    // 读取像片内外方系数文件
    ifstream file;
    file.open(path_ex);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
    string line;
    vector<vector<double>> ex; // 存储读取的数据
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
        ex.push_back(row);
    }
    file.close();

    photo_left = Photo(ex[0][4], ex[0][5], ex[0][6],ex[0][1]*pi/180,ex[0][2]*pi/180,ex[0][3]*pi/180);
    photo_right = Photo(ex[1][4], ex[1][5], ex[1][6],ex[1][1]*pi/180,ex[1][2]*pi/180,ex[1][3]*pi/180);
    
    file.open(path_in);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
    vector<vector<double>> in; // 存储读取的数据
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
        in.push_back(row);
    }
    file.close();
    photo_left.f = in[0][1]/1000;
    photo_left.x_0 = in[0][2]/1000;
    photo_left.y_0 = in[0][3]/1000;
    photo_left.m = in[0][4];

    photo_right.f = in[1][1]/1000;
    photo_right.x_0 = in[1][2]/1000;
    photo_right.y_0 = in[1][3]/1000;
    photo_right.m = in[1][4];

    vector<vector<double>> coordinate_l; // 存储读取的数据
    file.open(path_l);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
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
        coordinate_l.push_back(row);
    }
    file.close();

    vector<vector<double>> coordinate_r; // 存储读取的数据
    file.open(path_r);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
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
        coordinate_r.push_back(row);
    }

    for (int i = 0; i < coordinate_l.size(); i++)
    {
        Mat_<double> point(2, 1);
        point.at<double>(0, 0) = coordinate_l[i][4]/1000.0;
        point.at<double>(1, 0) = coordinate_l[i][3]/1000.0;
        photo_left.ctl_photo.addPoint(int(coordinate_l[i][0]), point);
    }

    for (int i = 0; i < coordinate_r.size(); i++)
    {
        Mat_<double> point(2, 1);
        point.at<double>(0, 0) = coordinate_r[i][4]/1000.0;
        point.at<double>(1, 0) = coordinate_r[i][3]/1000;
        photo_right.ctl_photo.addPoint(int(coordinate_r[i][0]), point);
    }
    return 0;
}


int datareader_3(string path_ex, string path_in, string path_l, string path_r, Photo& photo_left, Photo& photo_right)
{
    // 读取像片内外方系数文件
    ifstream file;
    file.open(path_ex);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
    string line;
    vector<vector<double>> ex; // 存储读取的数据
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
        ex.push_back(row);
    }
    file.close();

    photo_left = Photo(ex[0][4], ex[0][5], ex[0][6],ex[0][1]*pi/180,ex[0][2]*pi/180,ex[0][3]*pi/180);
    photo_right = Photo(ex[1][4], ex[1][5], ex[1][6],ex[1][1]*pi/180,ex[1][2]*pi/180,ex[1][3]*pi/180);
    
    file.open(path_in);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
    vector<vector<double>> in; // 存储读取的数据
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
        in.push_back(row);
    }
    file.close();
    photo_left.f = in[0][1]/1000;
    photo_left.x_0 = in[0][2]/1000;
    photo_left.y_0 = in[0][3]/1000;
    photo_left.m = in[0][4];

    photo_right.f = in[1][1]/1000;
    photo_right.x_0 = in[1][2]/1000;
    photo_right.y_0 = in[1][3]/1000;
    photo_right.m = in[1][4];

    vector<vector<double>> coordinate_l; // 存储读取的数据
    file.open(path_l);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
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
        coordinate_l.push_back(row);
    }
    file.close();

    vector<vector<double>> coordinate_r; // 存储读取的数据
    file.open(path_r);
    if (!file.is_open())
    {
        cout<< "fail to open file";
        return -1;
    }
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
        coordinate_r.push_back(row);
    }

    for (int i = 0; i < coordinate_l.size(); i++)
    {
        Mat_<double> point(2, 1);
        point.at<double>(0, 0) = coordinate_l[i][3]/1000.0;
        point.at<double>(1, 0) = coordinate_l[i][4]/1000.0;
        photo_left.ctl_photo.addPoint(int(coordinate_l[i][0]), point);
    }

    for (int i = 0; i < coordinate_r.size(); i++)
    {
        Mat_<double> point(2, 1);
        point.at<double>(0, 0) = coordinate_r[i][3]/1000.0;
        point.at<double>(1, 0) = coordinate_r[i][4]/1000;
        photo_right.ctl_photo.addPoint(int(coordinate_r[i][0]), point);
    }
    return 0;
}

map<int, GCP> readin_GCP(string file_path)
{
    ifstream gcp_file(file_path);
    if(!gcp_file.is_open())
    {
        cout << "Error: cannot open file " << file_path << endl;
    }
    // 读取gcp
    map<int, GCP> gcp_map;
    while(!gcp_file.eof())
    {
        GCP gcp;
        int id;
        int flag;
        double X, Y, Z;
        gcp_file >> id >> X >> Y >> Z >> flag;
        gcp.X = -Y;
        gcp.Y = -Z;
        gcp.Z = -X;
        gcp_map[id] = gcp;
    }
    return gcp_map;
}
// 读入observed_point
map<int, GCP> readin_op(string file_path, map<int, GCP> gcp_map)
{
    ifstream op_file(file_path);
    if(!op_file.is_open())
    {
        cout << "Error: cannot open file " << file_path << endl;
    }

    map<int, GCP> op_map;
    while(!op_file.eof())
    {
        GCP gcp;
        int id;
        double x, y;
        op_file >> id >> x >> y;
        x = 2128 - x;
        y = y - 1416;
        // x = x;
        // y = 4256 - y;
        if(gcp_map.find(id) == gcp_map.end())
        {
            continue;
        }
        gcp.x = x;
        gcp.y = y;
        gcp.X = gcp_map[id].X;
        gcp.Y = gcp_map[id].Y;
        gcp.Z = gcp_map[id].Z;
        op_map[id] = gcp;
    }
    return op_map;
}

map<int, UNP> readin_un(string file_path_l, string file_path_r)
{
    ifstream un_file_l(file_path_l);
    ifstream un_file_r(file_path_r);
    if(!un_file_l.is_open() || !un_file_r.is_open())
    {
        cout << "Error: cannot open file " << file_path_l << " or " << file_path_r << endl;
    }

    map<int, UNP> un_map;
    while(!un_file_l.eof() && !un_file_r.eof())
    {
        UNP un;
        int id;
        double x_l, y_l, x_r, y_r;
        un_file_l >> id >> x_l >> y_l;
        un_file_r >> id >> x_r >> y_r;

        un.x_l = 2128 -x_l;
        un.y_l = y_l- 1416;
        un.x_r = 2128 -x_r;
        un.y_r = y_r- 1416;
        un_map[id] = un;
    }
    return un_map;
}
