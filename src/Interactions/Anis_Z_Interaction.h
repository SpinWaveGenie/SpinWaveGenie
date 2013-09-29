#ifndef __AnisZ_Interaction_H__
#define __AnisZ_Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "Cell.h"
#include "Interaction.h"

class Anis_Z_Interaction: public Interaction
{
public:
    Anis_Z_Interaction(double value_in, std::string sl_r_in);
    void Update_Interaction(double value_in, std::string sl_r_in);
    void calcConstantValues(boost::shared_ptr<Cell> cell);
    void calcChangingValues(boost::shared_ptr<Cell> cell, Eigen::Vector3d K);
    void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN, int quadrant);
    std::vector<std::string> sublattices() const;
private:
    std::string sl_r;
    double value;
    int r,M;
    std::complex<double> LNrr,LNrrM,LNrMr,LNrMrM;
};

#endif
