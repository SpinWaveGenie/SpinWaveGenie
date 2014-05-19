#ifndef __ExchangeInteractionSameSublattice_H__
#define __ExchangeInteractionSameSublattice_H__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Cell/Cell.h"
#include "Interactions/Interaction.h"
#include "Containers/Matrices.h"
#include "Cell/Neighbors.h"


class ExchangeInteractionSameSublattice: public Interaction
{
public:
    ExchangeInteractionSameSublattice(std::string name, double value, std::string sl_r, double min, double max);
    void Update_Interaction(double value, std::string sl_r, double min, double max);
    virtual void updateValue(double value_in);
    virtual std::string getName();
    void calcConstantValues(Cell& cell);
    void calculateEnergy(Cell& cell, double &energy);
    void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~ExchangeInteractionSameSublattice(){};
private:
    Neighbors neighbors;
    std::string name,sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    std::complex<double> gamma_rs;
    std::complex<double> LNr,LNs,LNrr,LNss,LNrs,LNrsM;
};

#endif
