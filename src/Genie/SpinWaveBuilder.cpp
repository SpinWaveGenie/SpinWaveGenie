#include "SpinWaveBuilder.h"
#include <iostream>

using namespace std;

namespace SpinWaveGenie
{

    SpinWaveBuilder::SpinWaveBuilder(){}

    SpinWaveBuilder::SpinWaveBuilder(Cell& cellIn)
    {
        cell = cellIn;
    }
    
    void SpinWaveBuilder::updateCell(Cell & cellIn)
    {
        cell = cellIn;
    }

void SpinWaveBuilder::addInteraction(std::unique_ptr<Interaction> in)
{
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    interactions.push_back(in.release());
    interactions.sort();
}

void SpinWaveBuilder::addInteraction(Interaction* in)
{
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    //cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
    interactions.push_back(in);
    interactions.sort();
}
    
    void SpinWaveBuilder::updateInteraction(string name, double value)
    {
        for (auto it = interactions.begin();it!=interactions.end();it++)
        {
            if (name.compare(it->getName()) == 0)
                it->updateValue(value);
        }
    }

double SpinWaveBuilder::getEnergy()
{
    double energy = 0.0;
    for(auto it=interactions.begin();it!=interactions.end();it++)
    {
        it->calculateEnergy(cell,energy);
        //cout << it->getName() << " " << energy << endl;
    }
    //cout << endl;
    return energy;
}

Eigen::VectorXcd SpinWaveBuilder::getFirstOrderTerms()
{
    Eigen::VectorXcd firstOrder;
    firstOrder.setZero(2*cell.size());
    for (auto iter = interactions.begin(); iter != interactions.end(); iter++)
    {
        // vector<string> sls = iter->sublattices();
        // for(vector<string>::iterator iter2 = sls.begin();iter2 !=sls.end();++iter2)
        // {
        // cout << (*iter2) << " ";
        // }
        // cout << endl;
        //firstOrder.setZero(2*M);
        iter->calculateFirstOrderTerms(this->cell,firstOrder);
        //cout << firstOrder[2] << " " << firstOrder[8] << endl;
        //cout << firstOrder.transpose() << endl;
    }
    return firstOrder;
}

SpinWave SpinWaveBuilder::Create_Element()
{
    //cout << "cell check(Create_Element): " << cell.begin()->getName() << endl;
    for( auto in = interactions.begin();in!=interactions.end();in++)
    {
        in->calcConstantValues(cell);
    }
    SpinWave SW(cell,interactions);
    Eigen::VectorXcd firstOrder = getFirstOrderTerms();
    if (firstOrder.norm() > 0.01)
    {
        cout << "Warning! Nonzero first order terms present." << endl;
        cout << firstOrder.transpose() << endl;

    }
    return SW;
}
    
}
