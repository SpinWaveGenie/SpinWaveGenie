#ifndef __DMY_Interaction_H__
#define __DMY_Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Cell.h"
#include "Interaction.h"
#include "Cell/Neighbors.h"


class DM_Y_Interaction: public Interaction
{
public:
    DM_Y_Interaction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void Update_Interaction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void calcConstantValues(Cell& cell);
    void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN,int quadrant);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~DM_Y_Interaction(){};
private:
    Neighbors neighbors;
    std::string sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    double value0,value1,value2,value3;
    double z_rs;
    std::complex<double> gamma_rs;
};

#endif