// MDLW.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <conio.h>
#include <cstdlib>
//Argon
const long double sigma = 3.4 * powl(10, -10);
const long double m0 = 39.95 * 1.6747 * powl(10, -27);
const long double eps = 120 * 8.314;
//Cube in sigma
const int Lx = 1;
const int Ly = 1;
const int Lz = 1;
//
const int seed = 2021;
const double dt = 0.0001;
const double t0 = 0.0;
const double t1 = 0.0;
const int N = 2*2*2;
const int qN = 2;
const int vmax = 2;

//
const std::string dataFileName = "data.txt";    
void generateCubeUniformDistibutedPoins(
    float Lx, float Ly, float Lz, size_t qN, std::vector<double>& r0)
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
    float vv0;
    //1
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {

        }
        vv0 = -vmax + 2 * rand() %vmax;
    }
    //sum(p)->0

}
// Пока оставить
void save2file(const std::vector<double>& r0, size_t N)
{
    std::ofstream data;
    data.open(dataFileName,std::ios_base::out);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            data << r0[3 * i + j] << " ";
        }
        data << "\n";
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
    std::vector<double> r0;
    r0.resize(3*N);
    std::vector<double> v0;
    v0.resize(3 * N);
    generateCubeUniformDistibutedPoins(Lx, Ly, Lz, qN, r0);
    generateSpeed(vmax, v0);
    save2file(r0,N);
    plotScatter3d();
   
}