#include "dlt.h"

// l系数近似值解算
VectorXd my_DLT::init_l(vector<GCP> vector_gcp)
{
    VectorXd L_init;
	// 获取对应控制点
	
    MatrixXd A(ctl_point_index.size() * 2, 11);
    MatrixXd B(ctl_point_index.size() * 2, 1);
	A.fill(0);
    for (int i = 0; i < vector_gcp.size(); i++)
    {
        A(2 * i, 0) = vector_gcp[i].X;
        A(2 * i, 1) = vector_gcp[i].Y;
        A(2 * i, 2) = vector_gcp[i].Z;
        A(2 * i, 3) = 1;

        A(2 * i, 8) = vector_gcp[i].X * vector_gcp[i].x;
        A(2 * i, 9) = vector_gcp[i].Y * vector_gcp[i].x;
        A(2 * i, 10) = vector_gcp[i].Z * vector_gcp[i].x;

        A(2 * i + 1, 4) = vector_gcp[i].X;
        A(2 * i + 1, 5) = vector_gcp[i].Y;
        A(2 * i + 1, 6) = vector_gcp[i].Z;
        A(2 * i + 1, 7) = 1;

        A(2 * i + 1, 8) = vector_gcp[i].X * vector_gcp[i].y;
        A(2 * i + 1, 9) = vector_gcp[i].Y * vector_gcp[i].y;
        A(2 * i + 1, 10) = vector_gcp[i].Z * vector_gcp[i].y;

        B(2 * i, 0) = -vector_gcp[i].x;
        B(2 * i + 1, 0) = -vector_gcp[i].y;
    }
	L_init = (A.transpose() * A).inverse() * A.transpose() * B;
    return L_init;
}

VectorXd my_DLT::precise_l(vector<GCP> vector_gcp, VectorXd L_init)
{
    VectorXd L_precise = L_init;
    double x0, y0;
    double f;

    MatrixXd M(vector_gcp.size() * 2, 16);
    MatrixXd W(vector_gcp.size() * 2, 1);
	M.fill(0);
    VectorXd X;
    VectorXd V;
    for (int turn = 0; turn < 20; turn++) 
    {
        for (int i = 0; i < vector_gcp.size(); i++)
        {
            double A = L_precise[8] * vector_gcp[i].X + L_precise[9] * vector_gcp[i].Y + L_precise[10] * vector_gcp[i].Z + 1;
            x0 = -(L_precise[0] * L_precise[8] + L_precise[1] * L_precise[9] + L_precise[2] * L_precise[10]) / (L_precise[8] * L_precise[8] + L_precise[9] * L_precise[9] + L_precise[10] * L_precise[10]);
            y0 = -(L_precise[4] * L_precise[8] + L_precise[5] * L_precise[9] + L_precise[6] * L_precise[10]) / (L_precise[8] * L_precise[8] + L_precise[9] * L_precise[9] + L_precise[10] * L_precise[10]);
            double m = (L_precise(0, 0) * L_precise(0, 0) + L_precise(1, 0) * L_precise(1, 0) + L_precise(2, 0) * L_precise(2, 0)) / (L_precise(8, 0) * L_precise(8, 0) + L_precise(9, 0) * L_precise(9, 0) + L_precise(10, 0) * L_precise(10, 0)) - x0 * x0;
            double n = (L_precise(4, 0) * L_precise(4, 0) + L_precise(5, 0) * L_precise(5, 0) + L_precise(6, 0) * L_precise(6, 0)) / (L_precise(8, 0) * L_precise(8, 0) + L_precise(9, 0) * L_precise(9, 0) + L_precise(10, 0) * L_precise(10, 0)) - y0 * y0;
            double s = (L_precise(0, 0) * L_precise(4, 0) + L_precise(1, 0) * L_precise(5, 0) + L_precise(2, 0) * L_precise(6, 0)) / (L_precise(8, 0) * L_precise(8, 0) + L_precise(9, 0) * L_precise(9, 0) + L_precise(10, 0) * L_precise(10, 0)) - x0 * y0;
            f = sqrt((m * n - s * s) / n);

            double r = sqrt(pow(vector_gcp[i].x - x0, 2) + pow(vector_gcp[i].y - y0, 2));

            M(2 * i, 0) = -vector_gcp[i].X / A;
            M(2 * i, 1) = -vector_gcp[i].Y / A;
            M(2 * i, 2) = -vector_gcp[i].Z / A;
            M(2 * i, 3) = -1.0 / A;
            M(2 * i, 8) = -vector_gcp[i].X * vector_gcp[i].x / A;
            M(2 * i, 9) = -vector_gcp[i].Y * vector_gcp[i].x / A;
            M(2 * i, 10) = -vector_gcp[i].Z * vector_gcp[i].x / A;

            M(2 * i + 1, 4) = -vector_gcp[i].X / A;
            M(2 * i + 1, 5) = -vector_gcp[i].Y / A;
            M(2 * i + 1, 6) = -vector_gcp[i].Z / A;
            M(2 * i + 1, 7) = -1.0 / A;
            M(2 * i + 1, 8) = -vector_gcp[i].X * vector_gcp[i].y / A;
            M(2 * i + 1, 9) = -vector_gcp[i].Y * vector_gcp[i].y / A;
            M(2 * i + 1, 10) = -vector_gcp[i].Z * vector_gcp[i].y / A;

            M(2 * i, 11) = -(vector_gcp[i].x - x0) * pow(r, 2);
            M(2 * i, 12) = -(vector_gcp[i].x - x0) * pow(r, 4);
            M(2 * i, 13) = -(vector_gcp[i].x - x0) * pow(r, 6);
            M(2 * i, 14) = -(2 * pow(vector_gcp[i].x - x0, 2) + pow(r, 2));
            M(2 * i, 15) = -(2 * (vector_gcp[i].x - x0) * (vector_gcp[i].y - y0));

            M(2 * i + 1, 11) = -(vector_gcp[i].y - y0) * pow(r, 2);
            M(2 * i + 1, 12) = -(vector_gcp[i].y - y0) * pow(r, 4);
            M(2 * i + 1, 13) = -(vector_gcp[i].y - y0) * pow(r, 6);
            M(2 * i + 1, 14) = -(2 * (vector_gcp[i].x - x0) * (vector_gcp[i].y - y0));
            M(2 * i + 1, 15) = -(2 * pow(vector_gcp[i].y - y0, 2) + pow(r, 2));

            W(2 * i, 0) = vector_gcp[i].x / A;
            W(2 * i + 1, 0) = vector_gcp[i].y / A;
        }

        X = (M.transpose() * M).inverse() * M.transpose() * W;
        L_precise = X.head(11);
    }

    V = M * X - W;
    double sigma = sqrt((V.dot(V)) / (V.size() - 16));
    MatrixXd Q = (M.transpose() * M).inverse();
    Q = Q * sigma * sigma;
    cout << "sigma:" << sigma << endl;
    for (int i = 0; i < 11; i++) {
        cout << "l" << setprecision(5) << i + 1 << "\t" << X[i] << "\t" << sqrt(Q(i, i)) << endl;
    }
    cout << "像点观测值残差：" << endl;
    for (int i = 0; i < vector_gcp.size(); i++) {
        cout << ctl_point_index[i] << "\t" << V[2 * i] << "\t" << V[2 * i + 1] << endl;
    }

    cout << X[11] << "\t" << X[12] << "\t" << X[13] << "\t" << X[14] << "\t" << X[15] << "\t" << endl;

    return X;
}

void my_DLT::refine_point(vector<GCP> &vector_gcp, VectorXd L_precise)
{
    //4.待定点像点误差改正
    double x0 = -(L_precise[0] * L_precise[8] + L_precise[1] * L_precise[9] + L_precise[2] * L_precise[10]) / (L_precise[8] * L_precise[8] + L_precise[9] * L_precise[9] + L_precise[10] * L_precise[10]);
    double y0 = -(L_precise[4] * L_precise[8] + L_precise[5] * L_precise[9] + L_precise[6] * L_precise[10]) / (L_precise[8] * L_precise[8] + L_precise[9] * L_precise[9] + L_precise[10] * L_precise[10]);
	for (int i = 0; i < vector_gcp.size(); i++)
	{
        double r = sqrt(pow((vector_gcp[i].x - x0), 2) + pow((vector_gcp[i].y - y0), 2));//计算像点像径
        double tx, ty;
        tx = vector_gcp[i].x + ((vector_gcp[i].x - x0) * (L_precise(11) * pow(r, 2) + L_precise(12, 0) * pow(r, 4) + L_precise(13, 0) * pow(r, 6)) + L_precise(14, 0) * (pow(r, 2) + pow((vector_gcp[i].x - x0), 2)) + 2 * L_precise(15, 0) * (vector_gcp[i].x - x0) * (vector_gcp[i].y - y0));
        ty = vector_gcp[i].y + ((vector_gcp[i].y - y0) * (L_precise(11) * pow(r, 2) + L_precise(12, 0) * pow(r, 4) + L_precise(13, 0) * pow(r, 6)) + L_precise(14, 0) * (pow(r, 2) + pow((vector_gcp[i].y - y0), 2)) + 2 * L_precise(15, 0) * (vector_gcp[i].x - x0) * (vector_gcp[i].y - y0));
        vector_gcp[i].x = tx;
        vector_gcp[i].y = ty;
	}
}

void my_DLT::init_u(vector<UNP> &un_points, VectorXd LL, VectorXd RL)
{
    // 5.待定点物方空间坐标近似值的解算
	for (int i = 0; i < un_points.size() ; i++) 
    {
		
        MatrixXd A(4, 3);
        MatrixXd L(4, 1);
        VectorXd X;
        //求系数矩阵
        A(0, 0) = LL(0, 0) + un_points[i].x_l * LL(8, 0);
        A(0, 1) = LL(1, 0) + un_points[i].x_l * LL(9, 0);
        A(0, 2) = LL(2, 0) + un_points[i].x_l * LL(10, 0);

        A(1, 0) = LL(4, 0) + un_points[i].y_l * LL(8, 0);
        A(1, 1) = LL(5, 0) + un_points[i].y_l * LL(9, 0);
        A(1, 2) = LL(6, 0) + un_points[i].y_l * LL(10, 0);

        A(2, 0) = RL(0, 0) + un_points[i].x_r * RL(8, 0);
        A(2, 1) = RL(1, 0) + un_points[i].x_r * RL(9, 0);
        A(2, 2) = RL(2, 0) + un_points[i].x_r * RL(10, 0);

        A(3, 0) = RL(4, 0) + un_points[i].y_r * RL(8, 0);
        A(3, 1) = RL(5, 0) + un_points[i].y_r * RL(9, 0);
        A(3, 2) = RL(6, 0) + un_points[i].y_r * RL(10, 0);
        //求常数项矩阵
        L(0, 0) = -(LL(3, 0) + un_points[i].x_l);
        L(1, 0) = -(LL(7, 0) + un_points[i].y_l);
        L(2, 0) = -(RL(3, 0) + un_points[i].x_r);
        L(3, 0) = -(RL(7, 0) + un_points[i].y_r);

        X = (A.transpose() * A).inverse() * A.transpose() * L;
        un_points[i].X = X(0);
        un_points[i].Y = X(1);
        un_points[i].Z = X(2);
        cout<< "X: " << un_points[i].X << " Y: " << un_points[i].Y << " Z: " << un_points[i].Z << endl;
       
        cout << "\t" << X(0) << "\t" << X(1) << "\t" << X(2) << endl;
        
	}
}

//6.待定点物方空间坐标精确值的解算
void my_DLT::precise_u(vector<UNP> &un_points, VectorXd LL, VectorXd RL)
{
	for (int i = 0; i < un_points.size(); i++) 
    {
        MatrixXd A(4, 3);
        MatrixXd L(4, 1);
        VectorXd X;
        
        for (int turn = 0; turn < 10; turn++) 
        {
            double LA = LL(8, 0) * un_points[i].X + LL(9, 0) * un_points[i].Y + LL(10, 0) * un_points[i].Z + 1;
            double RA = RL(8, 0) * un_points[i].X + RL(9, 0) * un_points[i].Y + RL(10, 0) * un_points[i].Z + 1;

            //计算系数矩阵
            A(0, 0) = (LL(0, 0) + un_points[i].x_l * LL(8, 0))/LA;
            A(0, 1) = (LL(1, 0) + un_points[i].x_l * LL(9, 0))/LA;
            A(0, 2) = (LL(2, 0) + un_points[i].x_l * LL(10, 0))/LA;

            A(1, 0) = (LL(4, 0) + un_points[i].y_l * LL(8, 0))/LA;
            A(1, 1) = (LL(5, 0) + un_points[i].y_l * LL(9, 0))/LA;
            A(1, 2) = (LL(6, 0) + un_points[i].y_l * LL(10, 0))/LA;

            A(2, 0) = (RL(0, 0) + un_points[i].x_r * RL(8, 0))/RA;
            A(2, 1) = (RL(1, 0) + un_points[i].x_r * RL(9, 0))/RA;
            A(2, 2) = (RL(2, 0) + un_points[i].x_r * RL(10, 0))/RA;

            A(3, 0) = (RL(4, 0) + un_points[i].y_r * RL(8, 0))/RA;
            A(3, 1) = (RL(5, 0) + un_points[i].y_r * RL(9, 0))/RA;
            A(3, 2) = (RL(6, 0) + un_points[i].y_r * RL(10, 0))/RA;
            //求常数项矩阵
            L(0, 0) = -((LL(3, 0) + un_points[i].x_l))/LA;
            L(1, 0) = -((LL(7, 0) + un_points[i].y_l))/LA;
            L(2, 0) = -((RL(3, 0) + un_points[i].x_r))/RA;
            L(3, 0) = -((RL(7, 0) + un_points[i].y_r))/RA;

            // S = (N.transpose() * N).inverse() * N.transpose() * Q;
            X = (A.transpose() * A).inverse() * A.transpose() * L;
            un_points[i].X = X(0);
            un_points[i].Y = X(1);
            un_points[i].Z = X(2);
            cout<< "X: " << un_points[i].X << " Y: " << un_points[i].Y << " Z: " << un_points[i].Z << endl;
        }
        cout << "\t" << X(0) << "\t" << X(1) << "\t" << X(2) << endl;
	}
}

void my_DLT::cal_para(VectorXd L, string which)
{
    	
    double x0 = -(L[0] * L[8] + L[1] * L[9] + L[2] * L[10]) / (L[8] * L[8] + L[9] * L[9] + L[10] * L[10]);
    double y0 = -(L[4] * L[8] + L[5] * L[9] + L[6] * L[10]) / (L[8] * L[8] + L[9] * L[9] + L[10] * L[10]);
    double m = (L(0, 0) * L(0, 0) + L(1, 0) * L(1, 0) + L(2, 0) * L(2, 0)) / (L(8, 0) * L(8, 0) + L(9, 0) * L(9, 0) + L(10, 0) * L(10, 0)) - x0 * x0;
    double n = (L(4, 0) * L(4, 0) + L(5, 0) * L(5, 0) + L(6, 0) * L(6, 0)) / (L(8, 0) * L(8, 0) + L(9, 0) * L(9, 0) + L(10, 0) * L(10, 0)) - y0 * y0;
    double s = (L(0, 0) * L(4, 0) + L(1, 0) * L(5, 0) + L(2, 0) * L(6, 0)) / (L(8, 0) * L(8, 0) + L(9, 0) * L(9, 0) + L(10, 0) * L(10, 0)) - x0 * y0;
    double f = sqrt((m * n - s * s) / n);

    MatrixXd t(3, 3);
    Vector3d l, X;
    double r3 = 0;
    double a3 = 0, b3 = 0, c3 = 0;
    double A = 0, B = 0, C = 0;
    double ds = 0, db = 0, b1 = 0, b2 = 0;

    l(0, 0) = -L(3, 0);
    l(1, 0) = -L(7, 0);
    l(2, 0) = -1.0;

    t(0, 0) = L(0, 0);
    t(0, 1) = L(1, 0);
    t(0, 2) = L(2, 0);

    t(1, 0) = L(4, 0);
    t(1, 1) = L(5, 0);
    t(1, 2) = L(6, 0);

    t(2, 0) = L(8, 0);
    t(2, 1) = L(9, 0);
    t(2, 2) = L(10, 0);

    //求出像片的方向余弦
    r3 = sqrt(L(8, 0) * L(8, 0) + L(9, 0) * L(9, 0) + L(10, 0) * L(10, 0));
    a3 = L(8, 0) / r3;
    b3 = L(9, 0) / r3;
    c3 = L(10, 0) / r3;

    //求出比例尺不一误差ds和坐标轴不垂直性误差dβ
    A = (L(0, 0) * L(0, 0) + L(1, 0) * L(1, 0) + L(2, 0) * L(2, 0)) / (r3 * r3) - x0 * x0;
    B = (L(4, 0) * L(4, 0) + L(5, 0) * L(5, 0) + L(6, 0) * L(6, 0)) / (r3 * r3) - y0 * y0;
    C = (L(0, 0) * L(4, 0) + L(1, 0) * L(5, 0) + L(2, 0) * L(6, 0)) / (r3 * r3) - x0 * y0;
    ds = sqrt(A / B) - 1;
    db = asin(sqrt((C * C) / (B * A)));

    //求出方向余弦b1,b2
    b1 = (L(1, 0) / r3 + b3 * x0 + (y0 * b3 + L(5, 0)) * (1 + ds) * sin(db)) / f;
    b2 = (L(5, 0) / r3 + b3 * y0) * (1 + ds) * cos(db) / f;


    //求出三个外方位角元素phi,omega,kappa
    double phi, omega, kappa;
    phi = atan(-a3 / c3);
    omega = asin(-b3);
    kappa = atan(b1 / b2);

    //计算法方程式
    X = (t.transpose() * t).inverse() * t.transpose() * l;

    //计算三个外方位直线元素Xs,Ys,Zs
    double Xs, Ys, Zs;
    Xs = X(0, 0);
    Ys = X(1, 0);
    Zs = X(2, 0);

    cout << which<< "in_para:" <<"f:" <<f<< "\t" << "x0:" << x0 << "\t" << "y0:" << y0 << "\t" << "ds:" << ds << "\t" << "db:" << db << endl;
    cout << which<< "out_para" <<"Xs" << Xs << "\t" << "Ys" << Ys << "\t" << "Zs" << Zs << endl;
    cout << which<< "out_para" <<"phi" << phi << "\t" << "omega" << omega << "\t" << "kappa" << kappa << endl;
    cout << which<< "change_para" << L(11, 0) << "\t" << L(12, 0) << "\t" << L(13, 0) << "\t" << L(14, 0) << "\t" << L(15, 0) << endl; 
}

void my_DLT::cal_dist_to_52(vector<UNP> vector_unp)
{
    double x_52, y_52, z_52;
    for (int i = 0; i < vector_unp.size(); i++)
    {
        if (un_point_index[i] == 52)
        {
            x_52 = vector_unp[i].X;
            y_52 = vector_unp[i].Y;
            z_52 = vector_unp[i].Z;
        }
    }

    ofstream out("../result/dist_to_52.txt");
    cout << "与52号点距离：" << endl;
    for (int i = 0; i < vector_unp.size(); i++)
    {
        if (un_point_index[i] != 52)
        {
            double distance = sqrt(pow(vector_unp[i].X - x_52, 2) + pow(vector_unp[i].Y - y_52, 2) + pow(vector_unp[i].Z - z_52, 2));
            cout << un_point_index[i] << "\t" << distance << endl;
            out << un_point_index[i] << "\t" << distance << endl;
        }
    }
    out.close();
}


void my_DLT::dlt()
{
    vector<GCP> vector_gcp_l;
    vector<GCP> vector_gcp_r;
    vector<UNP> vector_unp;
    vector<UNP> vector_check;

    for (const auto& index : ctl_point_index) 
    {
        vector_gcp_l.push_back(gcps_l[index]);
        vector_gcp_r.push_back(gcps_r[index]);
    }

    
    for(const auto& pair : check_point_index)
    {
        vector_check.push_back(check_p[pair]);
    }

    for(const auto& pair : un_p)
    {
        vector_unp.push_back(pair.second);
    }

    // l系数近似值解算
    VectorXd LL_i = init_l(vector_gcp_l);
    VectorXd RL_i = init_l(vector_gcp_r);

    //2.内方位元素解算x0,y0解算
	double Lx0 = -(LL_i[0] * LL_i[8] + LL_i[1] * LL_i[9] + LL_i[2] * LL_i[10]) / (LL_i[8] * LL_i[8] + LL_i[9] * LL_i[9] + LL_i[10] * LL_i[10]);
	double Ly0 = -(LL_i[4] * LL_i[8] + LL_i[5] * LL_i[9] + LL_i[6] * LL_i[10]) / (LL_i[8] * LL_i[8] + LL_i[9] * LL_i[9] + LL_i[10] * LL_i[10]);
	double Rx0 = -(RL_i[0] * RL_i[8] + RL_i[1] * RL_i[9] + RL_i[2] * RL_i[10]) / (RL_i[8] * RL_i[8] + RL_i[9] * RL_i[9] + RL_i[10] * RL_i[10]);
	double Ry0 = -(RL_i[4] * RL_i[8] + RL_i[5] * RL_i[9] + RL_i[6] * RL_i[10]) / (RL_i[8] * RL_i[8] + RL_i[9] * RL_i[9] + RL_i[10] * RL_i[10]);
    cout << "Lx0: " << Lx0 << " Ly0: " << Ly0 << " Rx0: " << Rx0 << " Ry0: " << Ry0 << endl;

    // 3. l系数精确值的解算
    VectorXd LL_p = precise_l(vector_gcp_l, LL_i);
    VectorXd RL_p = precise_l(vector_gcp_r, RL_i);

    // 4. 待定点像点误差改正
    refine_point(vector_gcp_l, LL_p);
    refine_point(vector_gcp_r, RL_p);

    // 5. 待定点物方空间坐标近似值的解算
    init_u(vector_unp, LL_p, RL_p);
    init_u(vector_check, LL_p, RL_p);

    // 6. 待定点物方空间坐标精确值的解算
    precise_u(vector_unp, LL_p, RL_p);
    precise_u(vector_check, LL_p, RL_p);
    ofstream out2("../result/unp_p.txt");
    for(int i = 0; i < vector_unp.size(); i++)
    {
        cout << "un_point: " << un_point_index[i] <<vector_unp[i].X << " " << vector_unp[i].Y << " " << vector_unp[i].Z << endl;
        out2 << un_point_index[i] << " " << vector_unp[i].X << " " << vector_unp[i].Y << " " << vector_unp[i].Z << endl;
    }
    out2.close();
    ofstream out1("../result/check_p.txt");
    for(int i = 0; i < vector_check.size(); i++)
    {
        cout << "check point: " << check_point_index[i] <<vector_check[i].X << " " << vector_check[i].Y << " " << vector_check[i].Z << endl;
        out1 << check_point_index[i] << " " << vector_check[i].X << " " << vector_check[i].Y << " " << vector_check[i].Z << endl;
    }
    out1.close();

    // 7.内外方位元素的解算
    cal_para(LL_p, "left");
    cal_para(RL_p, "right");

    // 8.到52号点的距离计算
    cal_dist_to_52(vector_unp);
    // ofstream out1("../result/unp_p.txt");
    // ofstream out1("../result/clp_p.txt");
    // for (int i = 0; i < vector_unp.size(); i++) 
    // {
    //     out1 << un_point_index[i] << " " << vector_unp[i].X << " " << vector_unp[i].Y << " " << vector_unp[i].Z << endl;
    //     out1 << common_point_index[i]<< " " <<vector_unp[i].X << " " << vector_unp[i].Y << " " << vector_unp[i].Z << endl;
    // }
    // out1.close();
}


