#include "SpinWaveGenie/Genie/SpinWaveBuilder.h"
#include <iostream>

using namespace std;

namespace SpinWaveGenie
{

SpinWaveBuilder::SpinWaveBuilder() {}

SpinWaveBuilder::SpinWaveBuilder(Cell &cellIn) : cell(cellIn)
{
  for (const auto &iter : interactions) // r
  {
    interactions.push_back(iter->clone());
  }
}

void SpinWaveBuilder::updateCell(Cell &cellIn) { cell = cellIn; }

struct lessThanUniquePtr
{
  bool operator()(const unique_ptr<Interaction> &a, const unique_ptr<Interaction> &b) { return *a < *b; }
};

void SpinWaveBuilder::addInteraction(std::unique_ptr<Interaction> in)
{
  // cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
  // cout << "cell check(Add_Interaction): " << cell.begin()->getName() << endl;
  interactions.push_back(std::move(in));
  std::sort(interactions.begin(), interactions.end(), lessThanUniquePtr());
}

void SpinWaveBuilder::updateInteraction(string name, double value)
{
  for (auto & elem : interactions)
  {
    if (name.compare(elem->getName()) == 0)
      elem->updateValue(value);
  }
}

double SpinWaveBuilder::getEnergy()
{
  double energy = 0.0;
  for (auto & elem : interactions)
  {
    // energy = 0.0;
    elem->calculateEnergy(cell, energy);
    cout << elem->getName() << " " << energy / 4.0 << endl;
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
    elem->calculateFirstOrderTerms(this->cell, firstOrder);
    // cout << firstOrder[2] << " " << firstOrder[8] << endl;
    // cout << firstOrder.transpose() << endl;
  }
  return firstOrder;
}

SpinWave SpinWaveBuilder::createElement()
{
  // cout << "cell check(Create_Element): " << cell.begin()->getName() << endl;

  std::vector<std::unique_ptr<Interaction>> interactions_copy;
  for (const auto &iter : interactions) // r
  {
    interactions_copy.push_back(iter->clone());
  }

  for (auto &elem : interactions_copy)
  {
    elem->calcConstantValues(cell);
  }

  SpinWave SW(cell, std::move(interactions_copy));
  Eigen::VectorXcd firstOrder = getFirstOrderTerms();
  if (firstOrder.norm() > 0.1)
  {
    cout << "Warning! Nonzero first order terms present." << endl;
    cout << firstOrder.transpose() << endl;
  }
  return SW;
}
}
