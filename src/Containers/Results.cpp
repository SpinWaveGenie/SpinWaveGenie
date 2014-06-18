//
//  Results.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 6/1/14.
//
//

#include "Results.h"
#include <cmath>
using std::vector;

bool evalues_equal(const Point& a, const Point& b)
{
    // remove eigenvalues that are equal
    double EPS = 1.0e-5;
    return abs(a.frequency-b.frequency) < EPS;
}

void Results::clear()
{
    results.clear();
}

void Results::insert(Point value)
{
    results.push_back(value);
}

void Results::sort()
{
    std::sort(results.begin(),results.end());
}


void Results::uniqueSolutions()
{
    double EPS = 1.0e-5;
    int VP_pos;
    vector<Point> VI_unique;
    // Find unique eigenvalues
    VI_unique = results;
    std::sort(VI_unique.begin(),VI_unique.end());
    VI_unique.erase(unique(VI_unique.begin(),VI_unique.end(),evalues_equal),VI_unique.end());
    
    int NU = (int)VI_unique.size();
    VI_unique.resize(NU);
    
    for (int i=0;i<NU;i++)
    {
        VI_unique[i].intensity = 0.0;
    }
    
    for (int i=0;i<results.size();i++)
    {
        VP_pos = NU; //set position to a nonsense value
        for (int j=0;j<NU;j++)
        {
            if (std::abs(results[i].frequency - VI_unique[j].frequency) < EPS)
            {
                VP_pos = j;
                VI_unique[j].intensity += results[i].intensity;
                break;
            }
        }
        if (VP_pos== NU)
            std::cout << "error finding unique value" << std::endl;
    }
    
    results = VI_unique;
}

void Results::significantSolutions()
{
    double ETS = 0.001;
    vector<Point> VI_signif;
    
    for (int k=0;k!=results.size();k++)
    {
        if (results[k].intensity > ETS )
        {
            VI_signif.push_back(results[k]);
        }
    }
    results = VI_signif;
}