#include "SpinWaveGenie/Genie/SpinWaveBuilder.h"
#include <iostream>

using namespace std;

namespace SpinWaveGenie
{

SpinWaveBuilder::SpinWaveBuilder() {}

SpinWaveBuilder::SpinWaveBuilder(Cell &cellIn) { cell = cellIn; }

void SpinWaveBuilder::updateCell(Cell &cellIn) { cell = cellIn; }

void SpinWaveBuilder::addInteraction(std::unique_ptr<Interaction> in)
{
  // cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
  // cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
  interactions.push_back(in.release());
  interactions.sort();
}

void SpinWaveBuilder::updateInteraction(string name, double value)
{
  for (auto & elem : interactions)
  {
    if (name.compare(elem.getName()) == 0)
      elem.updateValue(value);
  }
}

double SpinWaveBuilder::getEnergy()
{
  double energy = 0.0;
  for (auto & elem : interactions)
  {
    // energy = 0.0;
    elem.calculateEnergy(cell, energy);
    cout << elem.getName() << " " << energy / 4.0 << endl;
  }
  // cout << endl;
  return energy;
}

Eigen::VectorXcd SpinWaveBuilder::getFirstOrderTerms()
{
  Eigen::VectorXcd firstOrder;
  firstOrder.setZero(2 * cell.size());
  for (auto & elem : interactions)
  {
    /*cout << iter->getName() << " ";
    vector<string> sls = iter->sublattices();
    for(vector<string>::iterator iter2 = sls.begin();iter2 !=sls.end();++iter2)
     {
     cout << (*iter2) << " ";
     }
    cout << endl;
    int M = cell.size();
    firstOrder.setZero(2*M);
    */
    elem.calculateFirstOrderTerms(this->cell, firstOrder);
    // cout << firstOrder[2] << " " << firstOrder[8] << endl;
    // cout << firstOrder.transpose() << endl;
  }
  return firstOrder;
}

SpinWave SpinWaveBuilder::createElement()
{
  // cout << "cell check(Create_Element): " << cell.begin()->getName() << endl;
  for (auto & elem : interactions)
  {
    elem.calcConstantValues(cell);
  }
  SpinWave SW(cell, interactions);
  Eigen::VectorXcd firstOrder = getFirstOrderTerms();
  if (firstOrder.norm() > 0.1)
  {
    cout << "Warning! Nonzero first order terms present." << endl;
    cout << firstOrder.transpose() << endl;
  }
  return SW;
}
}
