//
//  main.cpp
//  Spin Wave Fit
//
//  Created by Hahn, Steven E. on 1/7/13.
//  Copyright (c) 2013 Oak Ridge National Laboratory. All rights reserved.
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/version.hpp>
#include "SpinWave.h"
#include "Initializer.h"
#include "progressbar.h"

using namespace Eigen;
using namespace std;
using namespace boost;

struct mc_params
{
    double x0,y0,z0;
    double E_min,E_max;
    int E_points;
    shared_ptr<SW_Builder> builder;
};

double f(double x,double y,double z)
{
    Vector4d big_F,little_f,eval_points,weights;
    big_F << 0.3972,0.6295,-0.0314,0.0044;
    little_f << 13.2442,4.9034,0.3496,0.0;

    double f_Q = 0.0;
    for(int k=0;k<4;k++)
    {
        double s = 2.0*M_PI/5.4*sqrt(pow(x,2) + pow(y,2) + 0.5*pow(z,2))/(4.0*M_PI);
        f_Q += big_F[k]*exp(-1.0*little_f[k]*pow(s,2));
    }
    return f_Q;
}

VectorXd g(void *params)
{
    mc_params *exp_data = reinterpret_cast<mc_params*>(params);
    double x = exp_data->x0;
    double z = exp_data->z0;
    double E_min = exp_data->E_min;
    double E_max = exp_data->E_max;
    int E_points = exp_data->E_points;
    
    VectorXd integral,eval_points,weights;
    eval_points.setZero(32);
    weights.setZero(32);
    integral.setZero(E_points);
    
    double a,b,c;
    a = 1109.0;
    b = 0.0;
    c = 0.48;
    double norm = sqrt((a*c-b*b)/(a*M_PI*M_PI));
    
    eval_points << -7.12581390983,-6.40949814928,-5.81222594946,-5.27555098664,-4.77716450334,-4.30554795347,-3.85375548542,-3.41716749282,-2.99249082501,-2.57724953773,-2.16949918361,-1.76765410946,-1.37037641095,-0.97650046359,-0.584978765436,-0.194840741569,0.194840741569,0.584978765436,0.97650046359,1.37037641095,1.76765410946,2.16949918361,2.57724953773,2.99249082501,3.41716749282,3.85375548542,4.30554795347,4.77716450334,5.27555098664,5.81222594946,6.40949814928,7.12581390983;
    weights << 7.31067642754e-23, 9.23173653482e-19, 1.19734401957e-15, 4.21501019491e-13, 5.93329148347e-11, 4.09883215841e-09, 1.5741677944e-07, 3.65058512533e-06, 5.41658405999e-05, 0.000536268365495, 0.00365489032677, 0.0175534288315, 0.0604581309559, 0.151269734077, 0.277458142303, 0.375238352593, 0.375238352593, 0.277458142303, 0.151269734077, 0.0604581309559, 0.0175534288315, 0.00365489032677, 0.000536268365495, 5.41658405999e-05, 3.65058512533e-06, 1.5741677944e-07, 4.09883215841e-09, 5.93329148347e-11, 4.21501019491e-13, 1.19734401957e-15, 9.23173653482e-19, 7.31067642754e-23;
   
    for(int i=0;i!=E_points;i++)
    {
        double energy = E_min + (E_max-E_min)*(double)i/(double)(E_points-1);
        //cout << energy << endl;
        for(int j=0;j!=weights.size();j++)
        {
            double y = exp_data->y0 + eval_points[j]/sqrt(a);
            //cout << x << '\t' << y << '\t' << z << endl;
            double f_Q = f(x,y,z);
            SpinWave test = exp_data->builder->Create_Element(x,y,z);
            test.Calc();
            vector<double> frequencies = test.Get_Frequencies();
            vector<double> intensities = test.Get_Intensities();
            for(size_t k=0;k!=frequencies.size();k++)
            {
                //cout << frequencies[k] << endl;
                integral[i] += weights[j]*f_Q*f_Q*intensities[k]*exp(-2.0*b*(frequencies[k]-energy))*exp(-c*pow((frequencies[k]-energy),2));
            }
            
        }
    }
    integral = integral*norm;
    return integral;
}

VectorXd h(void *params)
{
    mc_params *exp_data = reinterpret_cast<mc_params*>(params);
    int E_points = exp_data->E_points;
    VectorXd integral,integral_sq;
    
#if BOOST_VERSION / 100 % 1000 >= 47
    boost::random::mt19937 gen;
    boost::random::uniform_real_distribution<double> die_x(1.8,2.2);
    boost::random::uniform_real_distribution<double> die_z(-3.2,-2.8);
#else
    boost::mt19937 gen;
    boost::uniform_real<double> dist_x(1.8,2.2);
    boost::uniform_real<double> dist_z(-3.2,-2.8);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > die_x(gen, dist_x);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > die_z(gen, dist_z);
#endif
    integral.setZero(E_points);
    for(int j=0;j!= 51;j++)
    {
#if BOOST_VERSION / 100 % 1000 >= 47
        exp_data->x0 = die_x(gen);
        exp_data->z0 = die_z(gen);
#else
        exp_data->x0 = die_x();
        exp_data->z0 = die_z();
#endif
        //cout << x << " " << y << " " << z << endl;
        VectorXd result = g(params);
        integral += result;
        //integral_sq += result*result;
    }
    integral = integral/51.0;
    /*double sigma = 0.0;
    for(int i=0;i!=integral.size();i++)
    {
        integral[i] = integral[i]/100001.0;
        integral_sq[i] = integral_sq[i]/100001.0;
        sigma += sqrt((integral_sq[i] - integral[i]*integral[i])/100000.0/100000.0);
        
    }
    cout << "sigma = " << sigma << endl;;
    */
    //cout << A << endl;
    return integral;
}

class ThreadClass {
public:
    int npoints,nproc,Epoints;
    MatrixXd figure;
    ThreadClass(int n) // Constructor
    {
        Eigen::initParallel();
        nproc = n;
        npoints = 64;
        Epoints = 51;
        figure.setZero(Epoints,npoints);
    }
    ~ThreadClass()
    {
        
    }
    // Destructor
    void Run(int i)
    {
        Init four_sl;
        four_sl.read_input("4sl_cell.xml");
        mc_params data;
        data.E_min = 0.0;
        data.E_max = 80.0;
        data.E_points = Epoints;
        data.builder = four_sl.get_builder();
        double x0,y0,z0,x1,y1,z1;
        x0=2.0;x1=2.0;
        y0=-1.5;y1=1.5;
        z0=-3.0;z1=-3.0;
        for(int m=i;m<npoints;m=m+nproc)
        {
            //cout << n << endl;
            data.x0 = x0 + (x1-x0)*m/(npoints-1);
            data.y0 = y0 + (y1-y0)*m/(npoints-1);
            data.z0 = z0 + (z1-z0)*m/(npoints-1);
            cout << data.x0 << " " << data.y0 << " " << data.z0 << endl;
            figure.col(m) = h((void *)&data);
        }
    }
};

int main(int argc, const char * argv[])
{
    /*Init four_sl;

    four_sl.read_input("4sl_cell.xml");
    //shared_ptr<SW_Builder> asdf = four_sl.get_builder();
    mc_params data;
    data.E_min = 0.0;
    data.E_max = 80.0;
    data.E_points = 51;
    data.builder = four_sl.get_builder();
    
    vector<double> E_array(data.E_points);
    for(int i=0;i<data.E_points;i++)
    {
        E_array[i] = data.E_min + (data.E_max-data.E_min)*(double)i/(double)(data.E_points-1);
    }
    
    double x0,y0,z0,x1,y1,z1;
    x0=2.0;x1=2.0;
    y0=0.0;y1=1.5;
    z0=-3.0;z1=-3.0;
    

    MatrixXd figure(data.E_points,npoints);
    ProgressBar pbar(npoints);
    pbar.start();
    for(int n=0;n!=npoints;n++)
    {
        pbar.update(n);
        //cout << n << endl;
        data.x0 = x0 + (x1-x0)*n/(npoints-1);
        data.y0 = y0 + (y1-y0)*n/(npoints-1);
        data.z0 = z0 + (z1-z0)*n/(npoints-1);
        //cout << x << " " << y << " " << z << endl;
        figure.col(n) = h((void *)&data);
    }
    pbar.finish();
    */
    
    boost::thread_group g;
    int n_threads = 12;
    ThreadClass tc(n_threads);
    for (int i=0;i<n_threads;i++)
    {
        boost::thread *t = new boost::thread(&ThreadClass::Run, &tc,i);
        g.add_thread(t);
    }
    g.join_all();
    
    std::ofstream file("test.txt");
    if (file.is_open())
    {
        file << tc.figure << '\n';
    }

    /*for(int i=0;i!=E_array.size();i++)
    {
        cout << E_array[i] << '\t' << intensities[i] << endl;
    }*/
   
    return 0;
}
