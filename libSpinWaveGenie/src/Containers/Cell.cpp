#include "SpinWaveGenie/Containers/Cell.h"
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdexcept>

using std::pair;
using std::string;
using std::cout;
using std::endl;

namespace SpinWaveGenie
{

struct CompareSublatticeNames
{
    CompareSublatticeNames(const std::string &name) : name(name) {}
    bool operator() (const Sublattice &arg)
    {
        return name.compare(arg.getName()) == 0;
    }
    std::string name;
};

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
    
    //basisVectors << a/2.0,a*sqrt(3.0)/2.0,0.0,
    //                a/2.0,-1.0*a*sqrt(3.0)/2.0,0.0,
    //                0.0,0.0,c;
    
    reciprocalVectors = 2.0*M_PI*basisVectors.inverse().transpose();

    //cout << "recip vectors equal" <<reciprocalVectors << endl;

}

void Cell::setBasisVectors(double scale, Matrix3 basis)
{
    basisVectors = scale*basis;
}

const Matrix3& Cell::getBasisVectors() const
{
    return basisVectors;
}

const Matrix3& Cell::getReciprocalVectors() const
{
    return reciprocalVectors;
}

void Cell::addSublattice(Sublattice& sl)
{
    std::string name = sl.getName();
    auto it = std::find_if(sublatticeInfo.begin(), sublatticeInfo.end(), CompareSublatticeNames(name));
    if (it != sublatticeInfo.end())
    {
        throw std::invalid_argument("sublattice already defined");
    }
    else
    {
        sublatticeInfo.push_back(sl);
    }
}

Sublattice& Cell::getSublattice(string name)
{
    auto it = std::find_if(sublatticeInfo.begin(), sublatticeInfo.end(), CompareSublatticeNames(name));
    if (it == sublatticeInfo.end())
    {
         throw std::invalid_argument("sublattice not found");
    }
    return *it;
}

    Sublattice& Cell::operator[](std::size_t position)
    {
        return sublatticeInfo[position];
    }

    
    
std::size_t Cell::getPosition(std::string name)
{
    auto it = std::find_if(sublatticeInfo.begin(), sublatticeInfo.end(), CompareSublatticeNames(name));
    if (it == sublatticeInfo.end())
    {
        throw std::invalid_argument("sublattice not found");
    }
    return std::distance(sublatticeInfo.begin(),it);
}

void Cell::addAtom(std::string name, double x, double y, double z)
{
    Vector3 scaled_position(x,y,z);
    
    //cout << "scaled= " << scaled_position.transpose() << endl;
    //cout << basisVectors << endl;
    
    Vector3 pos = scaled_position.transpose()*basisVectors;
    
    //cout << "unscaled= " << pos.transpose() << endl;
    //cout << " " <<endl;
    
    getSublattice(name).addAtom(pos[0], pos[1], pos[2]);
}

size_t Cell::size() const
{
    return sublatticeInfo.size();
}

Cell::Iterator Cell::begin()
{
    return sublatticeInfo.begin();
}

Cell::Iterator Cell::end()
{
    return sublatticeInfo.end();
}

Cell::ConstIterator Cell::cbegin()
{
    return sublatticeInfo.cbegin();
}

Cell::ConstIterator Cell::cend()
{
    return sublatticeInfo.cend();
}

}

