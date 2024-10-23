#include "space_resection.h"


void SpaceResection::resection()
{
    // 1. 构建A矩阵
    MatrixXd A(2*gcps.size(), 14);
    VectorXd L(2*gcps.size(), 1);
    MatrixXd bar(3, 1);
    MatrixXd R(3, 3);
    MatrixXd V;

    vector<GCP> vec_gcp;
    vector<int> vec_id;
    for (const auto& pair : gcps) 
    {
        vec_gcp.push_back(pair.second);
        vec_id.push_back(pair.first);
    }
    
    for(int times= 0; times<10; times++)
    {
        R(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
        R(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
        R(0, 2) = -sin(phi) * cos(omega);
        R(1, 0) = cos(omega) * sin(kappa);
        R(1, 1) = cos(omega) * cos(kappa);
        R(1, 2) = -sin(omega);
        R(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
        R(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
        R(2, 2) = cos(phi) * cos(omega);

        for (int i = 0; i < vec_gcp.size(); i++)
        {
            bar(0, 0) = vec_gcp[i].X - Xs;
            bar(1, 0) = vec_gcp[i].Y - Ys;
            bar(2, 0) = vec_gcp[i].Z - Zs;

            bar = R.inverse() * bar;

            A(2 * i + 0, 0) = 1 / bar(2, 0) * (R(0, 0) * f + R(0, 2) * (vec_gcp[i].x - x0));
            A(2 * i + 0, 1) = 1 / bar(2, 0) * (R(1, 0) * f + R(1, 2) * (vec_gcp[i].x - x0));
            A(2 * i + 0, 2) = 1 / bar(2, 0) * (R(2, 0) * f + R(2, 2) * (vec_gcp[i].x - x0));
            A(2 * i + 1, 0) = 1 / bar(2, 0) * (R(0, 1) * f + R(0, 2) * (vec_gcp[i].y - y0));
            A(2 * i + 1, 1) = 1 / bar(2, 0) * (R(1, 1) * f + R(1, 2) * (vec_gcp[i].y - y0));
            A(2 * i + 1, 2) = 1 / bar(2, 0) * (R(2, 1) * f + R(2, 2) * (vec_gcp[i].y - y0));
            A(2 * i + 0, 3) = (vec_gcp[i].y - y0) * sin(omega) - ((vec_gcp[i].x - x0) / f *
                ((vec_gcp[i].x - x0) * cos(kappa) - (vec_gcp[i].y - y0) * sin(kappa))
                + f * cos(kappa)) * cos(omega);
            A(2 * i + 0, 4) = -f * sin(kappa) - (vec_gcp[i].x - x0) / f *
                ((vec_gcp[i].x - x0) * sin(kappa) + (vec_gcp[i].y - y0) * cos(kappa));
            A(2 * i + 0, 5) = vec_gcp[i].y - y0;

            A(2 * i + 1, 3) = -(vec_gcp[i].x - x0) * sin(omega) - ((vec_gcp[i].y - y0) / f *
                ((vec_gcp[i].x - x0) * cos(kappa) - (vec_gcp[i].y - y0) * sin(kappa))
                - f * sin(kappa)) * cos(omega);
            A(2 * i + 1, 4) = -f * cos(kappa) - (vec_gcp[i].y - y0) / f *
                ((vec_gcp[i].x - x0) * sin(kappa) + (vec_gcp[i].y - y0) * cos(kappa));
            A(2 * i + 1, 5) = -(vec_gcp[i].x - x0);

            A(2 * i + 0, 6) = (vec_gcp[i].x - x0) / f;
            A(2 * i + 0, 7) = 1;
            A(2 * i + 0, 8) = 0;
            A(2 * i + 1, 6) = (vec_gcp[i].y - y0) / f;
            A(2 * i + 1, 7) = 0;
            A(2 * i + 1, 8) = 1;

            double r = sqrt(pow(vec_gcp[i].x - x0, 2) + pow(vec_gcp[i].y - y0, 2));

            double dx = (vec_gcp[i].x - x0) * (k1 * pow(r, 2) + k2 * pow(r, 4) + k3 * pow(r, 6))
                + p1 * (pow(r, 2) + 2 * pow(vec_gcp[i].x - x0, 2))
                + 2 * p2 * (vec_gcp[i].x - x0) * (vec_gcp[i].y - y0);
            double dy = (vec_gcp[i].y - y0) * (k1 * pow(r, 2) + k2 * pow(r, 4) + k3 * pow(r, 6))
                + p2 * (pow(r, 2) + 2 * pow(vec_gcp[i].y - y0, 2))
                + 2 * p1 * (vec_gcp[i].x - x0) * (vec_gcp[i].y - y0);

            A(2 * i + 0, 9) = -(vec_gcp[i].x - x0) * pow(r, 2);
            A(2 * i + 0, 10) = -(vec_gcp[i].x - x0) * pow(r, 4);
            A(2 * i + 0, 11) = -(vec_gcp[i].x - x0) * pow(r, 6);
            A(2 * i + 0, 12) = -(pow(r, 2) + 2 * pow(vec_gcp[i].x - x0, 2));
            A(2 * i + 0, 13) = -2 * (vec_gcp[i].x - x0) * (vec_gcp[i].y - y0);

            A(2 * i + 1, 9) = -(vec_gcp[i].y - y0) * pow(r, 2);
            A(2 * i + 1, 10) = -(vec_gcp[i].y - y0) * pow(r, 4);
            A(2 * i + 1, 11) = -(vec_gcp[i].y - y0) * pow(r, 6);
            A(2 * i + 1, 12) = -2 * (vec_gcp[i].x - x0) * (vec_gcp[i].y - y0);
            A(2 * i + 1, 13) = -(pow(r, 2) + 2 * pow(vec_gcp[i].y - y0, 2));

            L(2 * i + 0, 0) = vec_gcp[i].x + dx - x0 + f * bar(0, 0) / bar(2, 0);
            L(2 * i + 1, 0) = vec_gcp[i].y + dy - y0 + f * bar(1, 0) / bar(2, 0);
		}
        MatrixXd delta = (A.transpose() * A).inverse() * A.transpose() * L;
        cout<< delta << endl;
        Xs += delta(0, 0);
        Ys += delta(1, 0);
        Zs += delta(2, 0);
        phi += delta(3, 0);
        omega += delta(4, 0);
        kappa += delta(5, 0);

        f += delta(6, 0);
        x0 += delta(7, 0);
        y0 += delta(8, 0);

        k1 += delta(9, 0);
        k2 += delta(10, 0);
        k3 += delta(11, 0);
        p1 += delta(12, 0);
        p2 += delta(13, 0);
    }
    MatrixXd delta = (A.transpose() * A).inverse() * A.transpose() * L;
    V = A * delta - L;
    for(int i = 0; i < vec_gcp.size(); i++)
    {
        cout << "V" << vec_id[i] << ": " << V(2 * i + 0, 0) << " " << V(2 * i + 1, 0) << endl;
    }

    MatrixXd Q = (A.transpose() * A).inverse();
    double sigma = sqrt((V.transpose() * V)(0, 0) / (gcps.size() - 14));
    cout << "sigma:" << sigma << endl;
    cout << sigma*sqrt(Q(0,0))<< " " << sigma*sqrt(Q(1,1))<< " " << sigma*sqrt(Q(2,2))<< " " << sigma*sqrt(Q(3,3))<< " " << sigma*sqrt(Q(4,4))<< " " << sigma*sqrt(Q(5,5))<< " " << sigma*sqrt(Q(6,6))<< " " << sigma*sqrt(Q(7,7))<< " " << sigma*sqrt(Q(8,8))<< " " << sigma*sqrt(Q(9,9))<< " " << sigma*sqrt(Q(10,10))<< " " << sigma*sqrt(Q(11,11))<< " " << sigma*sqrt(Q(12,12))<< " " << sigma*sqrt(Q(13,13))<<endl;
    cout <<endl;
    
    cout << "Xs: " << Xs << endl;
    cout << "Ys: " << Ys << endl;
    cout << "Zs: " << Zs << endl;
    cout << "phi: " << phi << endl;
    cout << "omega: " << omega << endl;
    cout << "kappa: " << kappa << endl;
    cout << "f: " << f << endl;
    cout << "x0: " << x0 << endl;
    cout << "y0: " << y0 << endl;
    cout << "k1: " << k1 << endl;
    cout << "k2: " << k2 << endl;
    cout << "k3: " << k3 << endl;
    cout << "p1: " << p1 << endl;
    cout << "p2: " << p2 << endl;
    cout << endl; 
}