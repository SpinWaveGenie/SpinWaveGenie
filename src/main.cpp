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
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <thread>
#include "SpinWave.h"
#include "Initializer.h"
#include "progressbar.h"
#include "Cell/Neighbors.h"
#include "Cell/AtomIterator.h"
#include "extern/cubature.h"

using namespace Eigen;
using namespace boost;
namespace po = boost::program_options;
using namespace std;

struct mc_params
{
    double x0,y0,z0;
    double E_min,E_max;
    int E_points;
    SW_Builder builder;
};

int e_test(unsigned dim, const double *x, void *data,
           unsigned fdim, double *retval)
{
    mc_params *exp_data = reinterpret_cast<mc_params*>(data);

    double E_min = exp_data->E_min;
    double E_max = exp_data->E_max;
    int E_points = exp_data->E_points;

    Eigen::VectorXd integral;
    integral.setZero(E_points);
    
    double a,b,c;
    a = 1109.0;
    b = 0.0;
    c = 0.48;
    
    //double norm = sqrt((a*c-b*b)/(a*M_PI*M_PI));
    //double sigma = 1.0/sqrt(2.0*a);
    double sigma_energy = 1.0/sqrt(2.0*c);

    //sum += 1.0/exp(-a*pow(u,2));
    //cout << exp_data->x0 << '\t' << x[0] << '\t' << exp_data->y0 << '\t' << exp_data->z0 << endl;
    SpinWave test = exp_data->builder.Create_Element(exp_data->x0,x[0],exp_data->z0);
    test.Calc();
    vector<double> frequencies = test.Get_Frequencies();
    vector<double> intensities = test.Get_Intensities();
    
    double u = x[0] - exp_data->y0;
    for(size_t k=0;k!=frequencies.size();k++)
    {
        int min_bin = (int) (frequencies[k]-7.0*sigma_energy - E_min)*(E_points-1)/(E_max-E_min);
        int max_bin = (int) (frequencies[k]+7.0*sigma_energy - E_min)*(E_points-1)/(E_max-E_min);
        
        if (min_bin < 0)
            min_bin = 0;
        if (max_bin > E_points)
            max_bin = E_points;
        
        //for(int i=0;i!=E_points;i++)
        for(int i=min_bin;i!=max_bin;i++)
        {
            double energy = E_min + (E_max-E_min)*(double)i/(double)(E_points-1);
            
            //cout << energy << endl;
            //cout << frequencies[k] << endl;
            integral[i] += intensities[k]*exp(-c*pow(frequencies[k]-energy,2))*exp(-2.0*b*(frequencies[k]-energy)*u)*exp(-a*pow(u,2));
        }
    }
    
    for(int i=0;i!=E_points;i++)
    {
        retval[i] = integral[i];
    }
    return 0;
}

int f_test(unsigned dim, const double *x, void *data,
           unsigned fdim, double *retval)
{
    mc_params *exp_data = reinterpret_cast<mc_params*>(data);
    double y = exp_data->y0;
    //cout << x[0] << " " << x[1] << endl;
    exp_data->x0 = x[0];
    exp_data->z0 = x[1];
    //cout << exp_data->x0 << " " << exp_data->z0 << endl;
    int E_points = exp_data->E_points;
    
    Eigen::VectorXd integral;
    integral.setZero(E_points);
    
    double xmin, xmax;
    double tol, *err;

    tol = 1.0e-4;
    xmin = y - 0.2;
    xmax = y + 0.2;
    
    err = (double *) malloc(sizeof(double) * E_points);
    
    hcubature(E_points, e_test, (void *)exp_data,
              1, &xmin, &xmax,
              0, tol, 0, ERROR_INDIVIDUAL, retval, err);
    
    free(err);
    
    return 0;
}

class ThreadClass {
public:
    int npoints,nproc,Epoints;
    int pointsDone;
    MatrixXd figure;
    boost::mutex io_mutex; // The iostreams are not guaranteed to be thread-safe!
    ThreadClass(int n) // Constructor
    {
#if EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 1
        Eigen::initParallel();
#endif
        nproc = n;
        npoints = 101;
        pointsDone = 0;
        Epoints = 81;
        figure.setZero(Epoints,npoints);
    }
    ~ThreadClass()
    {
        
    }
    // Destructor
    void Run(int i, Init& four_sl)
    {
        boost::unique_lock<boost::mutex> scoped_lock(io_mutex);
        //Init four_sl(filename);
        //four_sl.read_input(filename);
        scoped_lock.unlock();
        mc_params data;
        data.E_min = 0.0;
        data.E_max = 80.0;
        data.E_points = Epoints;
        data.builder = four_sl.get_builder();
        
        double tol;
        unsigned dim, maxEval;
        
        dim = 2;
        tol = 1.0e-4;
        maxEval = 0;
        
        double x0,y0,z0,x1,y1,z1;
        x0=2.0;x1=2.0;
        y0=-1.5;y1=1.5;
        z0=-3.0;z1=-3.0;
        /*for(int m=i;m<npoints;m=m+nproc)
        {
            //cout << n << endl;
            data.x0 = x0 + (x1-x0)*m/(npoints-1);
            data.y0 = y0 + (y1-y0)*m/(npoints-1);
            data.z0 = z0 + (z1-z0)*m/(npoints-1);
            scoped_lock.lock();
            cout << data.x0 << " " << data.y0 << " " << data.z0 << endl;
            scoped_lock.unlock();
            figure.col(m) = h((void *)&data);
            //cout << figure.col(m) << endl;
        }*/
        
        vector<double> xmin(dim);
        vector<double> xmax(dim);
        
        xmin[0] = 1.8;
        xmax[0] = 2.2;
        xmin[1] = -3.2;
        xmax[1] = -2.8;
        
        for(int m=i;m<npoints;m=m+nproc)
        {
            data.x0 = x0 + (x1-x0)*m/(npoints-1);
            data.y0 = y0 + (y1-y0)*m/(npoints-1);
            data.z0 = z0 + (z1-z0)*m/(npoints-1);
            //cout << "Pos." << endl;
            //cout << x << " " << y << " " <<z << endl;
            //data.y0 = y;
            
            //scoped_lock.lock();
            //cout << data.x0 << " " << data.y0 << " " << data.z0 << endl;
            //scoped_lock.unlock();
            
            vector<double> val(Epoints);
            vector<double> err(Epoints);
        
            hcubature(Epoints, f_test, (void *)&data,
                      dim, &xmin[0], &xmax[0],
                      maxEval, tol, 0, ERROR_INDIVIDUAL, &val[0], &err[0]);
            
            //for(int i = 0;i!=Epoints;i++)
            //{
            //    double energy = Emin + (Emax-Emin)*(double)i/(double)(Epoints-1);
            //    cout << energy << "  " << val[i] << "  " << err[i] << endl;;
            //}
            
            for(int n=0;n<Epoints;n++)
            {
                figure(n,m) = val[n];
            }
            scoped_lock.lock();
            pointsDone++;
            scoped_lock.unlock();
            
            //cout << endl;
        }

    }
};

int main(int argc, char * argv[])
{
    string filename;
    int n_threads;
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("input", po::value<string>(), "set input filename")
        ("threads", po::value<int>(), "set number of threads")
        ;
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }
        if (vm.count("input")) {
            cout << "input filename was set to "
            << vm["input"].as<string>() << ".\n";
            filename = vm["input"].as<string>();
        } else {
            cout << "input filename was not set.\n";
            return 1;
        }
        if (vm.count("threads")) {
            cout << "Using "
            << vm["threads"].as<int>() << " processors.\n";
            n_threads = vm["threads"].as<int>();
        } else {
            cout << "number of threads was not set. Using 1 processor.\n";
            n_threads = 1;
            
        }
    }
    catch(std::exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
    
    Init four_sl(filename);
    
    boost::thread_group g;
    
    ThreadClass tc(n_threads);
    for (int i=0;i<n_threads;i++)
    {
        boost::thread *t = new boost::thread(&ThreadClass::Run, &tc, i, four_sl);
        g.add_thread(t);
    }
    
    ProgressBar pbar(tc.npoints);
    pbar.start();
    while(tc.pointsDone < tc.npoints)
    {
        sleep(1);
        pbar.update(tc.pointsDone);
    }
    pbar.finish();
    
    
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


