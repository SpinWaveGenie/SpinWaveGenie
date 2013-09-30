#include "Anis_Z_Interaction.h"

using namespace std;
using namespace Eigen;


Anis_Z_Interaction::Anis_Z_Interaction(double value_in, string sl_r_in)
{
    this->Update_Interaction(value_in, sl_r_in);
}

void Anis_Z_Interaction::Update_Interaction(double value_in, string sl_r_in)
{
    value = value_in;
    sl_r = sl_r_in;
}

vector<string> Anis_Z_Interaction::sublattices() const
{
    vector<string> sl;
    sl.push_back(sl_r);
    return sl;
}

void Anis_Z_Interaction::calcConstantValues(boost::shared_ptr<Cell> cell)
{
    //find location of r
    r= -1;
    M = 0;
    for (SublatticeIterator sl=cell->begin(); sl!=cell->end(); ++sl)
    {
        if ( sl_r == sl->getName())
            r = M;
        M++;
    }
    assert(r!=-1);

    double S = cell->getSublattice(sl_r).getMoment();
    double theta = cell->getSublattice(sl_r).getTheta();
    double X = value;
    
    LNrr = -0.5*X*S*(1.0-3.0*pow(cos(theta),2));
    LNrrM = -0.5*X*S*pow(sin(theta),2);
    LNrMr = -0.5*X*S*pow(sin(theta),2);
    LNrMrM = -0.5*X*S*(1.0-3.0*pow(cos(theta),2));
}

void Anis_Z_Interaction::calcChangingValues(boost::shared_ptr<Cell> cell, Vector3d K)
{
}

void Anis_Z_Interaction::Update_Matrix(Vector3d K, boost::shared_ptr<Cell> cell, MatrixXcd &LN, int quadrant)
{
    switch (quadrant)
    {
        case 0:
            LN(r,r) += LNrr;
            break;
        case 1:
            LN(r+M,r) += LNrMr;
            break;
        case 2:
            LN(r,r+M) += LNrrM;
            break;
        case 3:
            LN(r+M,r+M) += LNrMrM;
            break;
        default:
            //cout << "error: case must be between 0 and 3" << endl;
            break;
    }
}