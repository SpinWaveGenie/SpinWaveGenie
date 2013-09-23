#include "Cell.h"
#include "AtomIterator.h"
#include <boost/make_shared.hpp>

using namespace Eigen;
using namespace std;

bool Cell::FastCompare::operator ==(const FastCompare &other) const
{
    double error=1.0e-2;
    if (sl1 == other.sl1 && sl2 == other.sl2 && abs(min - other.min)<error && abs(max - other.max)<error)
        return true;
    else
        return false;
}

bool Cell::FastCompare::operator < ( const FastCompare &other) const
{
    double error=1.0e-2;
    bool answer;
    if (abs(max - other.max) < error)
    {
        if (abs(min - other.min) < error)
        {
            if (sl2 == other.sl2)
                answer = sl1->getName() < other.sl1->getName();
            else
                answer = sl2->getName() < other.sl2->getName();
        }
        else
            answer = min < other.min;
    }
    else
        answer = max < other.max;
    
    //std::cout << name1 << " " << name2 << " " << (name1 < name2) << answer << std::endl;
    
    return answer;
}

void Cell::set_basis_vectors(double a,double b, double c, double alpha_deg, double beta_deg, double gamma_deg)
{
    double alpha,beta,gamma;
    alpha = alpha_deg*M_PI/180.0;
    beta = beta_deg*M_PI/180.0;
    gamma = gamma_deg*M_PI/180.0;

    double ci,cj,ck;
    ci = c*cos(beta);
    cj = c*(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma);
    ck = c/sin(gamma)*sqrt(1.0-pow(cos(alpha),2)-pow(cos(beta),2)-pow(cos(gamma),2)+2.0*cos(alpha)*cos(beta)*cos(gamma));
    
    basis_vectors << a,0.0,0.0,
                     b*cos(gamma),b*sin(gamma),0.0,
                     ci,cj,ck;
    reciprocal_vectors = 2.0*M_PI*basis_vectors.inverse();

}

void Cell::set_basis_vectors(double scale, Eigen::Matrix3d basis)
{
    basis_vectors = scale*basis;
}

Eigen::Matrix3d Cell::get_basis_vectors()
{
    return basis_vectors;
}

Eigen::Matrix3d Cell::get_reciprocal_vectors()
{
    return reciprocal_vectors;
}

void Cell::add_sublattice(std::string name, boost::shared_ptr<Sublattice>& sl)
{
    sublattice_info.insert( pair<string,boost::shared_ptr<Sublattice> >(name,sl));
}

boost::shared_ptr<Sublattice> Cell::get_sublattice(string name)
{
    boost::shared_ptr<Sublattice> sublattice = sublattice_info[name];
    return sublattice;
}

void Cell::add_atom(std::string name, double x, double y, double z)
{
    Vector3d scaled_position;
    scaled_position << x,y,z;
    
    //cout << "scaled= " << scaled_position.transpose() << endl;
    //cout << basis_vectors << endl;
    
    Vector3d pos = basis_vectors*scaled_position;
    
    //cout << "unscaled= " << unscaled_position.transpose() << endl;
    
    sublattice_info[name]->addAtom(pos[0], pos[1], pos[2]);
}

vector<vector<double> >* Cell::get_neighbors(boost::shared_ptr<Sublattice>& sublattice1_in, boost::shared_ptr<Sublattice>& sublattice2_in , double min, double max)
{
    vector<vector<double> > results;
    Vector3d disp,disp_ang;
    FastCompare name;
    name.sl1 = sublattice1_in;
    name.sl2 = sublattice2_in;
    name.min = min;
    name.max = max;

    if (neighborCache.find(name) == neighborCache.end() )
    {
        //no benefit to iterating over the first sublattice. Hence we choose the first element
        AtomIterator atom1=sublattice1_in->begin();
        // Increase the size of the supercell until the list of neighbors does not change
        //   for two consecutive iterations. Number 5 suggested by Abinit code.
        //
        for (long supercellSize = 1;supercellSize<=5;supercellSize++)
        {
            //cout << supercellSize << endl;
            bool new_results = 0;
            {
            for (AtomIterator atom2=sublattice2_in->begin(); atom2!=sublattice2_in->end(); ++atom2)
            {
                for (long n1=-supercellSize;n1<=supercellSize;n1++)
                {
                    for (long n2=-supercellSize;n2<=supercellSize;n2++)
                    {
                        for (long n3=-supercellSize;n3<=supercellSize;n3++)
                        {
                            //find distance between supercells
                            disp << n1,n2,n3;
                            disp_ang = disp.transpose()*basis_vectors;
                            vector<double> temp(3);
                            temp[0] = (*atom1)[0] - (*atom2)[0] + disp_ang[0];
                            temp[1] = (*atom1)[1] - (*atom2)[1] + disp_ang[1];
                            temp[2] = (*atom1)[2] - (*atom2)[2] + disp_ang[2];
                            bool already_found = 0;
                            double norm = sqrt(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]);
                            if (norm < max && norm > min)
                            {
                                for (int k=0;k<results.size();k++)
                                {
                                    double difference = temp[0]-results[k][0] + temp[1]-results[k][1] + temp[2]-results[k][2];
                                    if (abs(difference) < 1.0e-5)
                                    {
                                        already_found = 1;
                                    }
                                }
                                if (already_found==0)
                                {
                                    new_results = 1;
                                    results.push_back(temp);
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
        neighborCache.insert(pair<FastCompare,vector<vector<double> > >(name,results) );
    }    
    //for (int i=0;i<neighborCache[name].size();i++)
    //{
    //    cout << neighborCache[name][i][0] << endl;
    //}
    //cout << "done" << endl;
    return &neighborCache[name];
}

int Cell::size()
{
    return sublattice_info.size();
}

SublatticeIterator Cell::begin()
{
    return SublatticeIterator(sublattice_info.begin());
}

SublatticeIterator Cell::end()
{
    return SublatticeIterator(sublattice_info.end());
}



