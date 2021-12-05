// MDLW.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <numeric>
#include <conio.h>
#include <cstdlib>
#include "Dense"
//alias
using Vector3d = Eigen::Vector3d;
//Argon
const long double sigma = 3.4 * powl(10, -10);
const long double m0 = 39.95 * 1.6747 * powl(10, -27);
const long double eps = 120 * 8.314;
//Cube in sigma
const Vector3d L{ 1, 1, 1 };
const std::vector<Vector3d> normal =
{
    Vector3d{0,0,0},
    Vector3d{0,0,1},
    Vector3d{0,0,-1},
    Vector3d{0,1,0},
    Vector3d{0,-1,0},
    Vector3d{1,0,0},
    Vector3d{-1,0,0}
};
const size_t NnormalCube = 6;
//
const int seed = 2021;
const double dt = 0.001;
const double t0 = 0.0;
const double t1 = 1.0;
const int N = 2*2*2;
const int qN = 2;
const int vmax = 2;
struct Entity
{
    Vector3d r;
    Vector3d v;
    Vector3d a;
};
//s
//
const std::string dataFileName = "data.txt";    
void generateCubeUniformDistibutedPoins(
    double Lx, double Ly, double Lz, size_t qN, std::vector<double>& r0)
{   
    
    double bx = Lx / (2 * qN);
    double by = Ly / (2 * qN);
    double bz = Lz / (2 * qN);
    size_t indx = 0;
    for (size_t i = 0; i < qN; i++)
    {
        for (size_t j = 0; j < qN; j++)
        {
            for (size_t k = 0; k < qN; k++)
            {
                indx = 3 * (qN * qN * i + qN * j + k);
                r0[indx] = bx + 2 * bx * i;
                r0[indx + 1] = by + 2 * by * j;
                r0[indx + 2] = bz + 2 * bz * k;
            }
        }
    }
}
void generateSpeed(int vmaxp, std::vector<double>& v0)
{
    std::array<double, 3> vsum{ 0 };
    std::array<double, 3> sump{ 0 };

    double r = 0;
    //1
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
            v0[3 * i + j] = -vmax + 2 * r* vmax;
            vsum[j] = vsum[j] + v0[3 * i + j];
        }

    }
    for (size_t j = 0; j < 3; j++)
    {
        vsum[j] = vsum[j] / N;
    }
    //+=or -=?
    //sum(p)->0 
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            v0[3 * i + j] = v0[3 * i + j] - vsum[j];
            sump[j] = sump[j] + v0[3 * i + j];
        }

    }
    std::cout << "Sum p =[";
    for (size_t j = 0; j < 3; j++)
    {
        std::cout << sump[j] << ",";
    }
    std::cout << "]\n";
}
void rij(double xi, double yi, double zi,
         double xj, double yj, double zj)
{
    ;
}
double rij(const Vector3d& ri, const Vector3d& rj)
{
    double max = sqrt(5);
    double dist;
    for (size_t i = 0; i < NnormalCube+1 ; i++)
    {
        //Достаточно квадрата? В книге посмотреть
        //dist = (ri - (rj + normal[i])).norm();
        dist = (ri - rj).norm();
        if (dist < max)
        {
            max = dist;
        }
    }
    return max;
}

Vector3d aij(const Vector3d& ri, const Vector3d& rj)
{
    double buff1 = 1 / rij(ri, rj);
    double buff2 = pow(buff1, 6.0);
    
    if (buff1 < 2.5)
    {
        return Vector3d{ 0,0,0 };
    }
    else
    {
        return -2 * (ri - rj) * buff1 * (2 * buff2 - 1) * buff2;
    }
}

void afull(const std::vector<double>& r, std::vector<double>& a)
{
    Vector3d aa{ 0,0,0 };
    // Создать отдельно или сам оптимизирует?
    Vector3d buff{ 0,0,0 };
    // Можно оптимизировать! до O(n^2/2)
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (i == j)
            {
                continue;
            }
            aa = aij(
                Vector3d{ r[3 * i],r[3 * i + 1],r[3 * i + 2] },
                Vector3d{ r[3 * j],r[3 * j + 1],r[3 * j + 2] }
            );
            buff = buff + aa;
        }
        a[3 * i] = buff[0];
        a[3 * i +1] = buff[1];
        a[3 * i + 2] = buff[2];
    }
}

void verleSpeed(std::vector<double>& r, std::vector<double>& v, std::vector<double>& a)
{
    auto _a = a;

    // 1 r(dr) -> r(dr+)
    for (size_t i = 0; i < 3 * N; i++)
    {
        r[i] = r[i] + v[i] * dt + 0.5 * a[i] * dt * dt;
    }
    // Create new f that involve in for 
    afull(r, a);
    for (size_t i = 0; i < 3 * N; i++)
    {
        v[i] = v[i] + 0.5 * (a[i] + _a[i]) * dt;
    }

}
struct Epoch
{
    double t = 0;
    std::vector<double> r;
    std::vector<double> v;
    std::vector<double> a;
};
void boundary1(std::vector<double>& r)
{
    // В книге взять
    // Проверить везде .x()!
    Vector3d buff{ 0, 0, 0 };
    //swith?
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            buff[j] = r[3 * i + j];
        }
        // 1 bottom
        if (buff[2] <0)
        {
            buff[2] = L.z();
        }
        // 2 up
        if (buff[2] > L.z())
        {
            buff[2] = 0;
        }
        // face
        if ((buff[2]>0)&&(buff[2] < L.z()))
        {
            if ((buff[0] > 0) && (buff[0] < L.x()))
            {
                // 3 west
                if (buff[1] < 0)
                {
                    buff[1] = L.y();

                }
                // 4 east
                if (buff[1] >L.y())
                {
                    buff[1] = 0;
                }
            }
            if ((buff[1] > 0) && (buff[1] < L.y()))
            {
                // 5 north
                if (buff[0] < 0)
                {
                    buff[0] = L.x();
                }
                // 6 south
                if (buff[0] > L.x())
                {
                    buff[0] = 0;
                }
            }

        }
        for (size_t j = 0; j < 3; j++)
        {
            r[3 * i + j] = buff[j];
        }
    }
}
std::vector<Epoch> task1(double t0, double t1, std::vector<double>& r0, std::vector<double>& v0)
{
    const size_t iters = static_cast<size_t>((t1 - t0) / dt);
    std::vector<double> a0;
    a0.resize(3 * N);
    // get a0
    afull(r0, a0);
    std::vector<Epoch> history;
    history.push_back(Epoch{ t0, r0, v0, a0 });
    for (size_t i = 1; i < iters; i++)
    {
        // Переименовать!
        verleSpeed(r0, v0, a0);
        history.push_back(Epoch{ t0 + i * dt, r0, v0, a0 });
        boundary1(r0);
    }
    return history;
}
double H(const std::vector<double>& r, const std::vector<double>& v)
{
    double T = 0;// Через среднюю скорость
    double U = 0;
    Vector3d vm{ 0,0,0 };
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            vm[j] = vm[j] + v[3 * i + j];
        }
    }
    T = (vm / N).norm();
    T = T * T;
    double buff1 = 0;
    double buff2 = 0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (i >=j)
            {
                continue;
            }
            buff1 = 1/ rij(
                Vector3d{ r[3 * i],r[3 * i + 1],r[3 * i + 2] },
                Vector3d{ r[3 * j],r[3 * j + 1],r[3 * j + 2] }
            );
            double ddd = rij(
                Vector3d{ r[3 * i],r[3 * i + 1],r[3 * i + 2] },
                Vector3d{ r[3 * j],r[3 * j + 1],r[3 * j + 2] }
            );
            if (buff1 < 2.5)
            {
                continue;
            }
            buff2 = pow(buff1, 6);
            U = U + (buff1 * buff1 - 1) * buff2;
        }

    }
    return T + U;
}
// Пока оставить
void save2file(const std::vector<double>& r0, size_t N, size_t dim = 3 )
{
    std::ofstream data;
    data.open(dataFileName,std::ios_base::out);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            data << r0[dim * i + j] << " ";
        }
        data << "\n";
    }
    data.close();
}       
void save2file(const std::vector<Epoch>& history)
{
    std::ofstream data;
    data.open("history.txt", std::ios_base::out);
    for (auto epoch : history)
    {
        ;
    }
    data.close();
}
void plotScatter2d()
{
    //system("gnuplot");
    std::string gnuplotPath("\"C:/Program Files/gnuplot/bin/gnuplot.exe\"");
    FILE* pipe = _popen(gnuplotPath.c_str(), "w");
    std::ifstream data;
    data.open("gnuplotCommandPlotData2d.txt", std::ios_base::in);
    std::string gnuplotCommand((std::istreambuf_iterator<char>(data)),
        std::istreambuf_iterator<char>());
    data.close();
    fprintf(pipe, gnuplotCommand.c_str());
    // Разобрать!
    fflush(pipe);
    _getch();
    _pclose(pipe);
}
void plotScatter3d()
{
    //system("gnuplot");
    std::string gnuplotPath("\"C:/Program Files/gnuplot/bin/gnuplot.exe\"");
    FILE* pipe = _popen(gnuplotPath.c_str(), "w");
    std::ifstream data;
    data.open("gnuplotCommandPlotData3d.txt", std::ios_base::in);
    std::string gnuplotCommand((std::istreambuf_iterator<char>(data)),
        std::istreambuf_iterator<char>());
    data.close();
    fprintf(pipe, gnuplotCommand.c_str());
    // Разобрать!
    fflush(pipe);
    _getch();
    _pclose(pipe);
}
int main()
{

    srand(seed);
    //Может массивы?
    std::cout << "initialization...\n";
    std::vector<double> r0;
    r0.resize(3*N);
    std::vector<double> v0;
    v0.resize(3 * N);
    generateCubeUniformDistibutedPoins(L.x(), L.y(), L.z(), qN, r0);
    generateSpeed(vmax, v0);
    std::cout << "calculate...\n";
    std::vector<Epoch> history = task1(t0, t1, r0, v0);
    std::cout << "calculate H...\n";
    std::vector<double> historyH;
    for (auto epoch : history)
    {
        historyH.push_back(H(epoch.r, epoch.v));
    }
    std::cout << "save... H\n";
    save2file(historyH, historyH.size(), 1);
    plotScatter2d();
}