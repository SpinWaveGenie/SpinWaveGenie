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
    // In principle, we only need to iterate over one atom in the first sublattice. However, iterating over
    // all atoms provides a good check that all atoms have the same number of neighbors in the same relative
    // positions
    for(Sublattice::Iterator atom1 = cell.getSublattice(sl1).begin(); atom1!=cell.getSublattice(sl1).end();++atom1)
    {
        // Increase the size of the supercell until the list of neighbors does not change
        // for two consecutive iterations. A 5x5x5 supercell should good enough for
        // any physical interaction. If not, a warning message will be printed.
        UniquePositions Neighbors;
        for (long supercellSize = 1;supercellSize<=5;supercellSize++)
        {
            //cout << supercellSize << endl;
            bool new_results = 0;
            for (Sublattice::Iterator atom2=cell.getSublattice(sl2).begin(); atom2!=cell.getSublattice(sl2).end(); ++atom2)
            {
                Vector3 atomdistance(atom2->get<0>()-atom1->get<0>(),atom2->get<1>()-atom1->get<1>(),atom2->get<2>()-atom1->get<2>());
                //cout << "atom2:  " << atom2->get<0>() << " " << atom2->get<1>() << " " << atom2->get<2>() << endl;
                //cout << "atom1:  " << atom1->get<0>() << " " << atom1->get<1>() << " " << atom1->get<2>() << endl;
                //cout << atomdistance.transpose() << endl;
                for (long n1=-supercellSize;n1<=supercellSize;n1++)
                {
                    for (long n2=-supercellSize;n2<=supercellSize;n2++)
                    {
                        for (long n3=-supercellSize;n3<=supercellSize;n3++)
                        {
                            //find distance between supercells
                            Vector3 dispRLU(n1,n2,n3);
                            Vector3 dispAng = dispRLU.transpose()*cell.getBasisVectors();
                            Vector3 relativeDistance = atomdistance + dispAng;
                            double norm = relativeDistance.norm();
                            // check if norm is between min and max
                            if (norm < max && norm > min)
                            {
                                Neighbors.insert(relativeDistance[0],relativeDistance[1],relativeDistance[2]);
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
        if (atom1 == cell.getSublattice(sl1).begin())
        {
            numberNeighbors = Neighbors.size();
            neighborList = Neighbors;
        }
        else
        {
            if (!(Neighbors == neighborList))
            {
                cout << "Old number of neighbors:" << numberNeighbors << endl;
                cout << "Old positions:" << endl;
                for(auto MyPosition = neighborList.begin(); MyPosition!= neighborList.end(); MyPosition++)
                {
                    cout << "\t" << MyPosition->get<0>() << " " << MyPosition->get<1>() << " " << MyPosition->get<2>() << endl;
                }
                cout << "New number of neighbors:" << Neighbors.size() << endl;
                for(auto MyPosition = Neighbors.begin(); MyPosition!= Neighbors.end(); MyPosition++)
                {
                    cout << "\t" << MyPosition->get<0>() << " " << MyPosition->get<1>() << " " << MyPosition->get<2>() << endl;
                }
            }
        }
    }
}

complex<double> Neighbors::getGamma(Vector3 K)
{
    complex<double> MXI (0.0,-1.0);
    complex<double> gamma_rs = complex<double>(0.0,0.0);
    for(Iterator nbr=this->begin(); nbr!=this->end(); ++nbr)
    {
        double dot_prod = K[0]*nbr->get<0>() + K[1]*nbr->get<1>() + K[2]*nbr->get<2>();
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



