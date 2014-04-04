#include <vector>
#include <iostream>
#include "IntegrateThetaPhi.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
int main()
{
    Matrix3d basisVectors;
    //basisVectors << 0.5,0.5,-0.5,
    //-0.5,0.5,0.5,
    //0.5,-0.5,0.5;
    
    //basisVectors << 0.5,0.5,0.0,
    //0.0,0.5,0.5,
    //0.5,0.0,0.5;
    
    //basisVectors << 3.0000000000,0.0000000000,0.0000000000,
    //-1.5000000000, 2.5980762114, 0.0000000000,
    //0.0000000000, 0.0000000000, 7.0000000000;
    
    basisVectors << 1.0,0.0,0.0,
    0.0,1.0,0.0,
    0.0,0.0,1.0;
    
    IntegrateThetaPhi func(0.0,5.0,401,basisVectors);
    
    for (int value = 0; value < 401; value++)
    {
        vector<double> result = func.getCut(0.0,0.0,(double)value/100.0);
        for ( auto it = result.begin(); it != result.end(); it++)
        {
            cout << (*it) << " ";
        }
        cout << endl;
    }
    return 0;
}
