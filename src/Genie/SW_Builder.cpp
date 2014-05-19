#include "SW_Builder.h"
#include <iostream>

using namespace std;

SW_Builder::SW_Builder()
{
}

SW_Builder::SW_Builder(Cell& cell_in)
{
    cell = cell_in;
    //cout << "cell check(Add_Cell): " << cell.begin()->getName() << endl;

}

void SW_Builder::addInteraction(std::unique_ptr<Interaction> in)
{
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    in->calcConstantValues(cell);
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    interactions.push_back(in.release());
    interactions.sort();
}

void SW_Builder::Add_Interaction(Interaction* in)
{
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    in->calcConstantValues(cell);
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    interactions.push_back(in);
    interactions.sort();
}

double SW_Builder::getEnergy()
{
    double energy = 0.0;
    for(auto it=interactions.begin();it!=interactions.end();it++)
    {
        it->calculateEnergy(cell,energy);
        cout << energy << endl;
    }
    return energy;
}

SpinWave SW_Builder::Create_Element()
{
    //cout << "cell check(Create_Element): " << cell.begin()->getName() << endl;
    SpinWave SW(cell,interactions);
    Eigen::VectorXcd firstOrder = SW.checkFirstOrderTerms();
    if (firstOrder.norm() > 0.01)
    {
        //cout << "Warning! Nonzero first order terms present." << endl;
        //cout << firstOrder.transpose() << endl;

    }
    return SW;
}
