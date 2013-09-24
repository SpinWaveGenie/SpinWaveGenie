#include "Anis_X_Interaction.h"

using namespace std;
using namespace Eigen;


Anis_X_Interaction::Anis_X_Interaction(double value_in, string sl_r_in)
{
    this->Update_Interaction(value_in, sl_r_in);
}

void Anis_X_Interaction::Update_Interaction(double value_in, string sl_r_in)
{
    value = value_in;
    sl_r = sl_r_in;
}

void Anis_X_Interaction::Update_Matrix(Vector3d K, boost::shared_ptr<Cell> cell, MatrixXcd &LN)
{
    complex<double> XI (0.0,1.0);
    //find location of r,s
    int r= -1;
    int M = 0;
    //CellIter sl(cell);
    //for(sl.First();!sl.IsDone();sl.Next())
    for (SublatticeIterator sl=cell->begin(); sl!=cell->end(); ++sl)
    {
        if ( sl_r == (*sl)->getName())
            r = M;
        M++;
    }
    assert(r!=-1);
        
    double S = cell->getSublattice(sl_r)->getMoment()[0];
    double theta = cell->getSublattice(sl_r)->getMoment()[1];
    double phi = cell->getSublattice(sl_r)->getMoment()[2];
    double X = value;
        
    LN(r,r)     -= 0.5*X*S*(pow(cos(theta),2)*pow(cos(phi),2)+pow(sin(phi),2)-2.0*pow(sin(theta),2)*pow(cos(phi),2));
    LN(r,r+M) -= 0.5*X*S*pow(cos(theta)*cos(phi)-XI*sin(phi),2);
    LN(r+M,r) -= 0.5*X*S*pow(cos(theta)*cos(phi)+XI*sin(phi),2);
    LN(r+M,r+M) -= 0.5*X*S*(pow(cos(theta),2)*pow(cos(phi),2)+pow(sin(phi),2)-2.0*pow(sin(theta),2)*pow(cos(phi),2));
}
