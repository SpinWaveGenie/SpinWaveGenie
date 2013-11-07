#include "SW_Builder.h"
#include <iostream>

using namespace std;


SW_Builder::SW_Builder()
{
}

SW_Builder::SW_Builder(Cell& cell_in)
{
    cell = cell_in;
}

void SW_Builder::Add_Interaction(Interaction* in)
{
    in->calcConstantValues(this->cell);
    interactions.push_back(in);
    interactions.sort();
}

Eigen::VectorXcd SW_Builder::checkFirstOrderTerms()
{
    int M=0;
    for (SublatticeIterator sl=cell.begin(); sl!=cell.end(); ++sl)
    {
        M++;
    }
    Eigen::VectorXcd firstOrder;
    firstOrder.setZero(2*M);
    return firstOrder;
}

SpinWave SW_Builder::Create_Element()
{
    SpinWave SW(cell,interactions);
    Eigen::VectorXcd firstOrder = SW.checkFirstOrderTerms();
    if (firstOrder.norm() > 0.001)
    {
        //cout << "Warning! Nonzero first order terms present." << endl;
    }
    return SW;
}
