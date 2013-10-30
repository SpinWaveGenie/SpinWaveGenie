#ifndef __Exch_Interaction_H__
#define __Exch_Interaction_H__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Cell.h"
#include "Interaction.h"

class Exch_Interaction: public Interaction
{
public:
    Exch_Interaction(double value, std::string sl_r,std::string sl_s, double min, double max);
    void Update_Interaction(double value, std::string sl_r,std::string sl_s, double min, double max);
    void calcConstantValues(Cell& cell);
    void calcChangingValues(Cell& cell, Eigen::Vector3d K);
    void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, Cell& cell, Eigen::MatrixXcd &LN, int quadrant);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~Exch_Interaction(){};
private:
    std::string sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    double Sr,Ss;
    Eigen::Matrix3d Frs,Fsr;
    double z_rs;
    std::complex<double> gamma_rs;
};

#endif
