#ifndef __ExchangeInteractionSameSublattice_H__
#define __ExchangeInteractionSameSublattice_H__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Genie/Neighbors.h"

namespace SpinWaveGenie
{

class ExchangeInteractionSameSublattice: public Interaction
{
public:
    ExchangeInteractionSameSublattice(std::string name, double value, std::string sl_r, double min, double max);
    void updateInteraction(double value, std::string sl_r, double min, double max);
    virtual void updateValue(double value_in);
    virtual const std::string& getName();
    void calcConstantValues(Cell& cell);
    void calculateEnergy(Cell& cell, double &energy);
    void calculateFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void updateMatrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~ExchangeInteractionSameSublattice(){};
private:
    Neighbors neighbors;
    std::string name,sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    std::complex<double> gamma_rs;
    std::complex<double> LNrr,LNss,LNrs,LNrsM;
};
}
#endif
