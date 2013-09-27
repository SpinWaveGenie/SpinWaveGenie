#include "Cell.h"
#include <iostream>
#include "AtomIterator.h"
#include <boost/make_shared.hpp>
#include <boost/functional/hash.hpp>

using namespace Eigen;
using namespace std;

bool Cell::NeighborsKey::operator ==(const NeighborsKey &other) const
{
    double error=1.0e-2;
    return (sl1.compare(other.sl1) == 0 && sl2.compare(other.sl2) == 0 && abs(min - other.min)<error && abs(max - other.max)<error);
}


std::size_t Cell::KeyHasher::operator()(const NeighborsKey& key) const
{
    std::size_t seed = 0;
    boost::hash_combine(seed, key.sl1);
    boost::hash_combine(seed, key.sl2);
    boost::hash_combine(seed, key.min);
    boost::hash_combine(seed, key.max);
    return seed;
}


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
    reciprocalVectors = 2.0*M_PI*basisVectors.inverse();

}

void Cell::setBasisVectors(double scale, Eigen::Matrix3d basis)
{
    basisVectors = scale*basis;
}

Eigen::Matrix3d Cell::getBasisVectors()
{
    return basisVectors;
}

Eigen::Matrix3d Cell::getReciprocalVectors()
{
    return reciprocalVectors;
}

void Cell::addSublattice(std::string name, boost::shared_ptr<Sublattice>& sl)
{
    sublatticeInfo.insert( pair<string,boost::shared_ptr<Sublattice> >(name,sl));
}

boost::shared_ptr<Sublattice> Cell::getSublattice(string name)
{
    boost::shared_ptr<Sublattice> sublattice = sublatticeInfo[name];
    return sublattice;
}

void Cell::addAtom(std::string name, double x, double y, double z)
{
    Vector3d scaled_position;
    scaled_position << x,y,z;
    
    //cout << "scaled= " << scaled_position.transpose() << endl;
    //cout << basisVectors << endl;
    
    Vector3d pos = basisVectors*scaled_position;
    
    //cout << "unscaled= " << unscaled_position.transpose() << endl;
    
    sublatticeInfo[name]->addAtom(pos[0], pos[1], pos[2]);
}

vector<vector<double> >* Cell::getNeighbors(string& sl1, string& sl2 , double min, double max)
{
    vector<vector<double> > results;
    Vector3d dispRLU,dispAng(3);
    NeighborsKey name;
    name.sl1 = sl1;
    name.sl2 = sl2;
    name.min = min;
    name.max = max;

    if (neighborCache.find(name) == neighborCache.end() )
    {
        //no benefit to iterating over the first sublattice. Hence we choose the first element
        AtomIterator atom1 = getSublattice(sl1)->begin();
        // Increase the size of the supercell until the list of neighbors does not change
        // for two consecutive iterations. A 5x5x5 supercell should good enough for
        // any physical interaction. if not a warning message will be printed.
        for (long supercellSize = 1;supercellSize<=5;supercellSize++)
        {
            //cout << supercellSize << endl;
            bool new_results = 0;
            {
            for (AtomIterator atom2=getSublattice(sl2)->begin(); atom2!=getSublattice(sl2)->end(); ++atom2)
            {
                for (long n1=-supercellSize;n1<=supercellSize;n1++)
                {
                    for (long n2=-supercellSize;n2<=supercellSize;n2++)
                    {
                        for (long n3=-supercellSize;n3<=supercellSize;n3++)
                        {
                            //find distance between supercells
                            dispRLU << n1,n2,n3;
                            dispAng = dispRLU.transpose()*basisVectors;
                            vector<double> relativeDistance(3);
                            relativeDistance[0] = (*atom2)[0] + dispAng[0] - (*atom1)[0];
                            relativeDistance[1] = (*atom2)[1] + dispAng[1] - (*atom1)[1];
                            relativeDistance[2] = (*atom2)[2] + dispAng[2] - (*atom1)[2];
                            bool already_found = 0;
                            double norm = sqrt(relativeDistance[0]*relativeDistance[0] + relativeDistance[1]*relativeDistance[1] + relativeDistance[2]*relativeDistance[2]);
                            // check if norm is between min and max
                            if (norm < max && norm > min)
                            {
                                // check if relativeDistance has already been calculated
                                for (int k=0;k<results.size();k++)
                                {
                                    double difference = relativeDistance[0]-results[k][0] + relativeDistance[1]-results[k][1] + relativeDistance[2]-results[k][2];
                                    if (abs(difference) < 1.0e-5)
                                    {
                                        already_found = 1;
                                        break;
                                    }
                                }
                                if (already_found==0)
                                {
                                    new_results = 1;
                                    results.push_back(relativeDistance);
                                }
                            }
                        }
                    }
                }
            }
            }
            if(!new_results)
                break;
            else if (supercellSize==5)
               cout << "Couldn't find all neighbors at specified distance" << endl;
        }
        neighborCache.insert(pair<NeighborsKey,vector<vector<double> > >(name,results) );
    }    
    //for (int i=0;i<neighborCache[name].size();i++)
    //{
    //    cout << neighborCache[name][i][0] << endl;
    //}
    //cout << "done" << endl;
    return &neighborCache[name];
}

size_t Cell::size()
{
    return sublatticeInfo.size();
}

SublatticeIterator Cell::begin()
{
    return SublatticeIterator(sublatticeInfo.begin());
}

SublatticeIterator Cell::end()
{
    return SublatticeIterator(sublatticeInfo.end());
}



