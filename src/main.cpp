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
#include <Eigen/Dense>
#include "SpinWavePlot.h"
#include "Initializer.h"
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>

using namespace boost;
namespace po = boost::program_options;
using namespace std;
using namespace Eigen;

class ThreadClass {
public:
    int npoints,nproc,points;
    double min,max;
    int pointsDone;
    MatrixXd figure;
    boost::mutex io_mutex; // The iostreams are not guaranteed to be thread-safe!
    ThreadClass(int n) // Constructor
    {
#if EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 1
        Eigen::initParallel();
#endif
        nproc = n;
        npoints = 401;
        pointsDone = 0;
        points = 321;
        figure.setZero(points,npoints);
        min = 0.0;
        max = 80.0;
    }
    ~ThreadClass()
    {
        
    }
    // Destructor
    void Run(int i, Init& four_sl)
    {
        boost::unique_lock<boost::mutex> scoped_lock(io_mutex);
        scoped_lock.unlock();
        SW_Builder builder = four_sl.get_builder();
        
        double tol;
        unsigned dim, maxEval;

        dim = 2;
        tol = 1.0e-4;
        maxEval = 0;
        
        double x0,y0,z0,x1,y1,z1;
        x0=2.0;x1=2.0;
        y0=-1.5;y1=1.5;
        z0=-3.0;z1=-3.0;
        
        vector<double> xmin(dim);
        vector<double> xmax(dim);
        
        xmin[0] = 1.8;
        xmax[0] = 2.2;
        xmin[1] = -3.2;
        xmax[1] = -2.8;
        
        IntegrateAxes tmp(builder,min,max,points);

        for(int m=i;m<npoints;m=m+nproc)
        {
            double x = x0 + (x1-x0)*m/(npoints-1);
            double y = y0 + (y1-y0)*m/(npoints-1);
            double z = z0 + (z1-z0)*m/(npoints-1);
            //cout << "Pos." << endl;
            //scoped_lock.lock();
            //cout << x << " " << y << " " << z << endl;
            //scoped_lock.unlock();
            vector<double> val = tmp.getCut(x,y,z);
            for(int n=0;n<points;n++)
            {
                figure(n,m) = val[n];
            }
            scoped_lock.lock();
            pointsDone++;
            scoped_lock.unlock();
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
    boost::unique_lock<boost::mutex> scoped_lock(tc.io_mutex);
    double npoints = tc.npoints;
    scoped_lock.unlock();
    boost::progress_display show_progress(npoints);
    int pointsDone = 0;
    while(pointsDone < npoints)
    {
        sleep(1);
        scoped_lock.lock();
        int diff = tc.pointsDone - pointsDone;
        pointsDone = tc.pointsDone;
        //cout << pointsDone << endl;
        scoped_lock.unlock();
        show_progress += diff;
    }
    g.join_all();
    
    std::ofstream file("test.txt");
    if (file.is_open())
    {
        file << tc.figure << '\n';
    }
    
    return 0;
}

