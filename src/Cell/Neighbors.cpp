#include <iostream>
#include "Neighbors.h"
#include "Cell.h"

using std::string;
using std::complex;
using std::cout;
using std::endl;

Neighbors::Neighbors()
{
    
}

void Neighbors::findNeighbors(Cell& cell, string& sl1, string& sl2 , double min, double max)
{
    
    Sublattice::Iterator atom1 = cell.getSublattice(sl1).begin();
    //no benefit to iterating over the first sublattice. Hence we choose the first element
    // Increase the size of the supercell until the list of neighbors does not change
    // for two consecutive iterations. A 5x5x5 supercell should good enough for
    // any physical interaction. if not a warning message will be printed.
    for (long supercellSize = 1;supercellSize<=5;supercellSize++)
    {
        //cout << supercellSize << endl;
        bool new_results = 0;
        for (Sublattice::Iterator atom2=cell.getSublattice(sl2).begin(); atom2!=cell.getSublattice(sl2).end(); ++atom2)
        {
            for (long n1=-supercellSize;n1<=supercellSize;n1++)
            {
                for (long n2=-supercellSize;n2<=supercellSize;n2++)
                {
                    for (long n3=-supercellSize;n3<=supercellSize;n3++)
                    {
                        //find distance between supercells
                        Vector3 dispRLU(n1,n2,n3);
                        Vector3 dispAng = dispRLU.transpose()*cell.getBasisVectors();
                        Vector3 relativeDistance = (*atom2) + dispAng - (*atom1);
                        bool already_found = 0;
                        double norm = relativeDistance.norm();
                        // check if norm is between min and max
                        if (norm < max && norm > min)
                        {
                            // check if relativeDistance has already been calculated
                            // cout << relativeDistance[0] << " " << relativeDistance[1] << " " << relativeDistance[2] << endl;
                            for (int k=0;k<neighborList.size();k++)
                            {
                                Vector3 difference = relativeDistance - neighborList[k];
                                if (difference.squaredNorm() < 1.0e-5)
                                {
                                    already_found = 1;
                                    break;
                                }
                            }
                            if (already_found==0)
                            {
                                new_results = 1;
                                neighborList.push_back(relativeDistance);
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
    
    numberNeighbors = neighborList.size();
    
    /*
     cout << name.sl1 << " " << name.sl2 << endl;
     for (int i=0;i<neighborCache[name].size();i++)
     {
     cout << neighborCache[name][i][0] << " " << neighborCache[name][i][1] << " " <<neighborCache[name][i][2] << endl;
     }
     cout << "done" << endl;
    */
}

complex<double> Neighbors::getGamma(Vector3 K)
{
    complex<double> MXI (0.0,-1.0);
    complex<double> gamma_rs = complex<double>(0.0,0.0);
    for(Iterator nbr=this->begin(); nbr!=this->end(); ++nbr)
    {
        double dot_prod = K.dot(*nbr);
        gamma_rs += exp(MXI*dot_prod);
    }
    return gamma_rs/numberNeighbors;
}

double Neighbors::getNumberNeighbors()
{
    return numberNeighbors;
}

Neighbors::Iterator Neighbors::begin()
{
    return Iterator(neighborList.begin());
}

Neighbors::Iterator Neighbors::end()
{
    return Iterator(neighborList.end());
}



