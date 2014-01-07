#include "Cell/Cell.h"
#include "Containers/Matrices.h"
#include <iostream>
#include <stdexcept>

using std::pair;
using std::string;
using std::cout;
using std::endl;

void Cell::setBasisVectors(double a,double b, double c, double alpha_deg, double beta_deg, double gamma_deg)
{
    //! <a href=https://github.com/mantidproject/documents/blob/master/Design/UBMatriximplementationnotes.pdf> Reference </a>
    double alpha,beta,gamma;
    alpha = alpha_deg*M_PI/180.0;
    beta = beta_deg*M_PI/180.0;
    gamma = gamma_deg*M_PI/180.0;

    double ci,cj,ck;
    ci = c*cos(beta);
    cj = c*(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma);
    ck = c/sin(gamma)*sqrt(1.0-pow(cos(alpha),2)-pow(cos(beta),2)-pow(cos(gamma),2)+2.0*cos(alpha)*cos(beta)*cos(gamma));
    
    basisVectors << a,0.0,0.0,
                    b*cos(gamma),b*sin(gamma),0.0,
                    ci,cj,ck;
    
    //cout << "basis vectors equal" <<basisVectors << endl;
    reciprocalVectors = 2.0*M_PI*basisVectors.inverse();
}

void Cell::setBasisVectors(double scale, Matrix3 basis)
{
    basisVectors = scale*basis;
}

const Matrix3 Cell::getBasisVectors()
{
    return basisVectors;
}

const Matrix3 Cell::getReciprocalVectors()
{
    return reciprocalVectors;
}

void Cell::addSublattice(Sublattice& sl)
{
    string name = sl.getName();
    sublatticeInfo.insert(pair<string,Sublattice>(name,sl));
}

Sublattice& Cell::getSublattice(string name)
{
    return sublatticeInfo[name];
}

const int Cell::getPosition(std::string name)
{
    int r = 0;
    for (Iterator sl=begin(); sl!=end(); ++sl)
    {
        if ( name.compare(sl->getName()) == 0)
            return r;
        r++;
    }
    throw std::invalid_argument("sublattice not found");
}

void Cell::addAtom(std::string name, double x, double y, double z)
{
    Vector3 scaled_position;
    scaled_position << x,y,z;
    
    //cout << "scaled= " << scaled_position.transpose() << endl;
    //cout << basisVectors << endl;
    
    Vector3 pos = basisVectors*scaled_position;
        
    sublatticeInfo[name].addAtom(pos[0], pos[1], pos[2]);
}

const size_t Cell::size()
{
    return sublatticeInfo.size();
}

Cell::Iterator Cell::begin()
{
    return Iterator(sublatticeInfo.begin());
}

Cell::Iterator Cell::end()
{
    return Iterator(sublatticeInfo.end());
}



