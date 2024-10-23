#ifndef DLT
#define DLT

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <random>
#include <chrono>

#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Core"

#include "data_reader.h"

using namespace std;
using namespace Eigen;


class my_DLT
{
public:
    map<int, GCP> gcps_l;
    map<int, GCP> gcps_r;
    map<int, UNP> check_p;
    map<int, UNP> un_p;
    vector<int> ctl_point_index;
    vector<int> check_point_index;
    vector<int> un_point_index;

    my_DLT(map<int, GCP> gcps_l, map<int, GCP> gcps_r, map<int, UNP> un_points)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);

        this->gcps_l = gcps_l;
        this->gcps_r = gcps_r;
        this->un_p = un_points;
        ofstream out("../result/check_o.txt");
        for (const auto& pair : gcps_l) 
        {
            if (gcps_r.find(pair.first) != gcps_r.end()) 
            {
                if(dis(gen) < 0.2)
                {
                    check_point_index.push_back(pair.first);
                    check_p[pair.first] = {gcps_l[pair.first].x, gcps_l[pair.first].y, gcps_r[pair.first].x, gcps_r[pair.first].y, 
                                            gcps_l[pair.first].X, gcps_l[pair.first].Y, gcps_l[pair.first].Z};
                    out << pair.first << " " << gcps_l[pair.first].X << " " << gcps_l[pair.first].Y << " " << gcps_l[pair.first].Z << endl;
                }
                else
                {
                    ctl_point_index.push_back(pair.first);
                }
            }
        }
        out.close();
        
        for (const auto& pair : un_points) 
        {
            un_point_index.push_back(pair.first);
        }
    }

    VectorXd init_l(vector<GCP> vector_gcp);
    VectorXd precise_l(vector<GCP> vector_gcp, VectorXd L_init);
    void refine_point(vector<GCP> &vector_gcp, VectorXd L_precise);
    void init_u(vector<UNP> &un_points, VectorXd LL, VectorXd RL);
    void precise_u(vector<UNP> &un_points, VectorXd LL, VectorXd RL);
    void cal_para(VectorXd L, string which);
    void cal_dist_to_52(vector<UNP> vector_unp);

    void dlt();
    



};

#endif