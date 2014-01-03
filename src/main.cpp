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
#include "SpinWave.h"
#include "Initializer.h"
#include "Cell/Neighbors.h"
#include "PointsAlongLine.h"
#include "Positions.h"

using namespace boost;
using namespace std;

int main(int argc, char * argv[])
{
    
    //Init four_sl;
    Init four_sl("/Users/svh/Documents/spin_wave_genie/examples/MnV2O4_cubic.xml");
    SW_Builder builder = four_sl.get_builder();
        
    //cout << check << endl;
    
    /*Cell cell = four_sl.get_cell();
    
    string sl_r = "Mn0";
    string sl_s = "Mn1";
    double min = 0.0;
    double max = 7.0;
    
    Neighbors neighborList;
    neighborList.findNeighbors(cell,sl_r,sl_s,min,max);
    size_t z_rs = neighborList.getNumberNeighbors();
    
    cout << "z_rs= " << z_rs << endl;
    for(Neighbors::Iterator nbr=neighborList.begin();nbr!=neighborList.end();++nbr)
    {
        cout << (*nbr)[0] << " " << (*nbr)[1] << " " << (*nbr)[2] << endl;
    }*/
    
    SpinWave test = builder.Create_Element();
    
    double SA = 2.3;
    double SB = 0.9;
    double theta = M_PI - 35.0*M_PI/180.0;
    double DB = -6.62711;
    double JBB = -9.80542;
    double DBz = 0.073016;
    
    for(int dbi = 0; dbi != 1; dbi++)
    {
    double JBBP = 6.56457*(1.0+(double)dbi/50.0);
    double JAB = SB*((6.0*JBB+6.0*JBBP+DB-3.0*DBz)*cos(theta)*sin(theta) -sqrt(2.0)*DB*(2.0*pow(cos(theta),2)-1.0))/(-9.0*SA*sin(theta));
    cout << "JBBP= " << JBBP << " " << " JAB= " << JAB << endl;
    test.updateValue("Jbbp",JBBP);
    test.updateValue("Jab",JAB);

    std::cout << "hello world" << std::endl;
    

    PointsAlongLine Line;
    Line.setFirstPoint(1.0,1.0,0.0);
    Line.setFinalPoint(3.0,3.0,0.0);
    Line.setNumberPoints(11);
    Positions KPoints = Line.getPoints();
        
    for(Positions::Iterator it = KPoints.begin(); it != KPoints.end(); it++)
    {
        double x = it->get<0>();
        double y = it->get<1>();
        double z = it->get<2>();

        //cout << "Pos." << endl;
        cout << x << " " << y << " " << z << " ";// << endl;
        test.createMatrix(x,y,z);
        test.Calc();
        vector<point> pts = test.getPoints();
        //cout << "Freq.  Int." << endl;
        for(vector<point>::iterator it2 = pts.begin();it2!=pts.end();it2++)
        {
            cout << (*it2).frequency << "  " << (*it2).intensity*10.0 << " " ;//<< endl;
        }
        cout << endl;
    }
    }
    return 0;
}
