#ifndef __DMZ_Interaction_H__
#define __DMZ_Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Cell.h"
#include "Interaction.h"

class DM_Z_Interaction: public Interaction
{
public:
    DM_Z_Interaction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void Update_Interaction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void calcConstantValues(Cell& cell);
    void calcChangingValues(Cell& cell, Eigen::Vector3d K);
    void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, Cell& cell, Eigen::MatrixXcd &LN,int quadrant);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~DM_Z_Interaction(){};
private:
    std::string sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    double tmp0,tmp1,tmp2,tmp3,tmp4;
    double z_rs;
    std::complex<double> gamma_rs;
};

#endif