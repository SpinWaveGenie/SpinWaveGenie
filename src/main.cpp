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
#include "SpinWave.h"
#include "Initializer.h"
#include "Cell/Neighbors.h"

using namespace boost;
using namespace std;

int main(int argc, char * argv[])
{
    //Init four_sl;
    Init four_sl("/Users/svh/Documents/spin_wave_genie/examples/YFeO3TwoSublattice.xml");
    SW_Builder builder = four_sl.get_builder();
        
    //cout << check << endl;
    
    Cell cell = four_sl.get_cell();
    
    string sl_r = "Spin0";
    string sl_s = "Spin1";
    double min = 0.0;
    double max = 1.0;
    
    Neighbors neighborList;
    neighborList.findNeighbors(cell,sl_r,sl_s,min,max);
    size_t z_rs = neighborList.getNumberNeighbors();
    
    cout << "z_rs= " << z_rs << endl;
    for(Neighbors::Iterator nbr=neighborList.begin();nbr!=neighborList.end();++nbr)
    {
        cout << (*nbr)[0] << " " << (*nbr)[1] << " " << (*nbr)[2] << endl;
    }

    int npoints = 11;
    double x,y,z,x0,y0,z0,x1,y1,z1;
    x0=0.0;x1=1.0;
    y0=0.0;y1=0.0;
    z0=0.0;z1=0.0;
    for(int m=0;m<npoints;m++)
    {
        //sleep(1);
        if(npoints==1)
        {
            x = x0 ;
            y = y0 ;
            z = z0 ;
        }
        else
        {
            x = x0 + (x1-x0)*m/(npoints-1);
            y = y0 + (y1-y0)*m/(npoints-1);
            z = z0 + (z1-z0)*m/(npoints-1);
        }
        //cout << "Pos." << endl;
        cout << x << " " << y << " " << z << " ";// << endl;
        SpinWave test = builder.Create_Element();
        test.createMatrix(x,y,z);
        test.Calc();
        vector<double> frequencies = test.Get_Frequencies();
        vector<double> intensities = test.Get_Intensities();
        //cout << "Freq.  Int." << endl;
        vector<double>::iterator it2 = intensities.begin();
        for(vector<double>::iterator it = frequencies.begin();it!=frequencies.end();it++)
        {
            cout << (*it) << "  ";// << (*it2) ;//<< endl;
            it2++;
        }
        cout << endl;
    }
    return 0;
}
