#ifndef __ExchangeInteraction_H__
#define __ExchangeInteraction_H__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "SpinWaveGenie/Cell/Cell.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Cell/Neighbors.h"

namespace SpinWaveGenie
{

class ExchangeInteraction: public Interaction
{
public:
    ExchangeInteraction(std::string name, double value, std::string sl_r,std::string sl_s, double min, double max);
    void updateInteraction(double value, std::string sl_r,std::string sl_s, double min, double max);
    virtual void updateValue(double value_in);
    virtual const std::string& getName();
    void calcConstantValues(Cell& cell);
    void calculateEnergy(Cell& cell, double &energy);
    void calculateFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void updateMatrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~ExchangeInteraction(){};
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
