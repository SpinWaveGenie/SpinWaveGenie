#ifndef __Exch_Interaction_H__
#define __Exch_Interaction_H__

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "Cell.h"
#include "Interaction.h"

class Exch_Interaction: public Interaction
{
public:
    Exch_Interaction(double value, std::string sl_r,std::string sl_s, double min, double max);
    void Update_Interaction(double value, std::string sl_r,std::string sl_s, double min, double max);
    void calcConstantValues(boost::shared_ptr<Cell> cell);
    void calcChangingValues(boost::shared_ptr<Cell> cell, Eigen::Vector3d K);
    void checkFirstOrderTerms(boost::shared_ptr<Cell> cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN, int quadrant);
    std::vector<std::string> sublattices() const;
private:
    std::string sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    double Sr,Ss,X;
    std::complex<double> G1,G2;
    Eigen::Matrix3d Frs;
    double z_rs;
    std::complex<double> gamma_rs;
};

#endif
