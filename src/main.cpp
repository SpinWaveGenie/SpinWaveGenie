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
#include "IntegrateThetaPhi.h"
#include "Initializer.h"
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>

using namespace boost;
namespace po = boost::program_options;
using namespace std;
using namespace Eigen;



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
    
    //cout << check << endl;
        
    Cell cell = four_sl.get_cell();
        
    string sl_r = "Mn0";
    string sl_s = "Mn";
    double min = 0.0;
    double max = 7.0;
        
        
        
    for (int i=0;i<2;i++)
    {
        string sl_si = sl_s + to_string(i);
        cout << sl_si << endl;
        Neighbors neighborList;
        neighborList.findNeighbors(cell,sl_r,sl_si,min,max);
        size_t z_rs = neighborList.getNumberNeighbors();
            
        cout << "z_rs= " << z_rs << endl;
        for(Neighbors::Iterator nbr=neighborList.begin();nbr!=neighborList.end();++nbr)
        {
            cout << sqrt(pow(nbr->get<0>(),2) + pow(nbr->get<1>(),2) +pow(nbr->get<2>(),2)) << " ";
            cout << nbr->get<0>() << " " << nbr->get<1>() << " " << nbr->get<2>() << endl;
        }
            
    }
    return 0;
}

