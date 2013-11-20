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
#include <boost/atomic.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>

using namespace boost;
namespace po = boost::program_options;
using namespace std;
using namespace Eigen;

class ThreadClass {
public:
    boost::atomic<int> npoints,nproc,points;
    boost::atomic<double> min,max;
    boost::atomic<int> pointsDone,nextPoint;
    MatrixXd figure;
    boost::mutex io_mutex; // The iostreams are not guaranteed to be thread-safe!
    ThreadClass(int n) // Constructor
    {
#if EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 1
        Eigen::initParallel();
#endif
        nproc.store(n);
        npoints.store(26);
        pointsDone.store(0);
        nextPoint.store(0);
        points.store(21);
        figure.setZero(points.load(),npoints.load());
        min.store(0.0);
        max.store(80.0);
    }
    ~ThreadClass()
    {
        
    }
    // Destructor
    void Run(int i, Init& four_sl)
    {
        boost::unique_lock<boost::mutex> scoped_lock(io_mutex);
        SW_Builder builder = four_sl.get_builder();
        scoped_lock.unlock();

        double x0,y0,z0,x1,y1,z1;
        x0=2.0;x1=2.0;
        y0=-1.5;y1=1.5;
        z0=-3.0;z1=-3.0;
        
        TwoDimGaussian resinfo;
        resinfo.a = 1109.0;
        resinfo.b = 0.0;
        resinfo.c = 0.48;
        resinfo.tol = 1.0e-2;
        resinfo.direction = 1;
        resinfo.SW = builder.Create_Element();
        
        axes_info axesinfo;
        axesinfo.x = true;
        axesinfo.y = true;
        axesinfo.z = true;
        axesinfo.dx = 0.2;
        axesinfo.dy = 0.05;
        axesinfo.dz = 0.2;
        axesinfo.tol = 1.0e-2;
        
        
        /*double x0,y0,z0,x1,y1,z1;
        x0=3.0;x1=3.0;
        y0=0.0;y1=0.0;
        z0=-4.0;z1=4.0;
         
        TwoDimGaussian resinfo;
        resinfo.a = 579.7;
        resinfo.b = -20.0;
        resinfo.c = 1.3;
        resinfo.tol = 1.0e-1;
        resinfo.direction = 2;
        resinfo.builder = builder;
         
        axes_info axesinfo;
        axesinfo.x = true;
        axesinfo.y = true;
        axesinfo.z = true;
        axesinfo.dx = 0.2;
        axesinfo.dy = 0.2;
        axesinfo.dz = 0.05;
        axesinfo.tol = 1.0e-1;
        */
        
        double tmin = min.load();
        double tmax = max.load();
        int tpoints = points.load();
        int tnpoints = npoints.load();

        TwoDimensionResolutionFunction res(resinfo, tmin,tmax, tpoints);
        IntegrateAxes tmp(axesinfo,res,tmin,tmax,tpoints);
        //for(int m=i;m<tnpoints;m=m+tnproc)
        while (true)
        {
            int m = nextPoint.load();

            if (m >= tnpoints)
            {
                break;
            }
            //cout << m << endl;
            nextPoint.fetch_add(1);
            double x = x0 + (x1-x0)*m/(tnpoints-1);
            double y = y0 + (y1-y0)*m/(tnpoints-1);
            double z = z0 + (z1-z0)*m/(tnpoints-1);
            //cout << "Pos." << endl;
            //scoped_lock.lock();
            //cout << "**** " << x << " " << y << " " << z << endl;
            //scoped_lock.unlock();
            vector<double> val = tmp.getCut(x,y,z);
            for(int n=0;n<tpoints;n++)
            {
                figure(n,m) = val[n];
            }
            pointsDone.fetch_add(1);
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
    scoped_lock.unlock();
    double npoints = tc.npoints;
    boost::progress_display show_progress(npoints);
    int pointsDone = 0;
    while(pointsDone < npoints)
    {
        sleep(1);
        int diff = tc.pointsDone.load() - pointsDone;
        if (diff > 0)
        {
            scoped_lock.lock();
            show_progress += diff;
            scoped_lock.unlock();
            pointsDone = tc.pointsDone.load();
        }
        //cout << pointsDone << endl;
    }
    g.join_all();
    
    std::ofstream file("test.txt");
    if (file.is_open())
    {
        file << tc.figure << '\n';
    }
    return 0;
}

