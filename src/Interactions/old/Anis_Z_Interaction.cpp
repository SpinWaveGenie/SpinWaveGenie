#include "Anis_Z_Interaction.h"

using namespace std;
using namespace Eigen;

Anis_Z_Interaction::Anis_Z_Interaction()
{
    
}

void Anis_Z_Interaction::Add_Interaction(double value_in, string sl_r_in)
{
    Anis_Z_Parameters in;
    in.value = value_in;
    in.sl_r = sl_r_in;
    exch_array.push_back(in);
}

void Anis_Z_Interaction::Update_Matrix(Vector3d K, boost::shared_ptr<Cell> cell, MatrixXcd &LN)
{
    
    std::vector<Anis_Z_Parameters>::iterator iter;
    for(iter=Anis_Z_array.begin();iter!=Anis_Z_array.end();iter++)
    {
        //find location of r,s
        int r= -1;
        
        CellIter sl(cell);
        for(sl.First();!sl.IsDone();sl.Next())
        {
            if ( iter->sl_r == sl.CurrentItem()->get_name())
                r = M;
            M++;
        }
        assert(r!=-1);
        
        double S = cell->get_sublattice(iter->sl_r)->get_moment()[0];
        double theta = cell->get_sublattice(iter->sl_r)->get_moment()[1];
        double X = iter->value;
               
        LN(r,r) -= 0.5*X*S*(1.0-3.0*pow(cos(theta),2));
        LN(r+M,r+M) -= 0.5*X*S*(1.0-3.0*pow(cos(theta),2));
        LN(r,r+M) -= 0.2*X*S*pow(sin(theta),2);
        LN(r+M,r) -= 0.2*X*S*pow(sin(theta),2);
    }
}