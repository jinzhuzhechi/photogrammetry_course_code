#ifndef SPACE_RESECTION
#define SPACE_RESECTION

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>

#include "data_reader.h"

using namespace std;
using namespace Eigen;



class SpaceResection
{
public:
    map<int, GCP> gcps;
    double f, x0, y0;
    double k1, k2, k3, p1, p2;
    double Xs, Ys, Zs;
    double phi, omega, kappa;

    void resection();

    SpaceResection(double f, double x_0, double y_0, 
                    double Xs, double Ys, double Zs, double phi, double omega, double kappa, 
                    double k1 = 0, double k2 = 0, double k3 = 0, double p1 = 0, double p2 = 0)
    {
        this->f = f;
        this->x0 = x_0;
        this->y0 = y_0;
        this->Xs = Xs;
        this->Ys = Ys;
        this->Zs = Zs;
        this->phi = phi;
        this->omega = omega;
        this->kappa = kappa;
        this->k1 = k1;
        this->k2 = k2;
        this->k3 = k3;
        this->p1 = p1;
        this->p2 = p2;
    }

    void addGCP(int id, double x, double y, double X, double Y, double Z)
    {
        GCP gcp;
        gcp.x = x;
        gcp.y = y;
        gcp.X = X;
        gcp.Y = Y;
        gcp.Z = Z;
        this->gcps[id] = gcp;
    }
};


#endif