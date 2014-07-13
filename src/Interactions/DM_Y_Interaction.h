#ifndef __DMY_Interaction_H__
#define __DMY_Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Cell/Cell.h"
#include "Interactions/Interaction.h"
#include "Cell/Neighbors.h"

namespace SpinWaveGenie
{
    
class DM_Y_Interaction: public Interaction
{
public:
    DM_Y_Interaction(std::string name, double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void updateInteraction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in);
    void updateValue(double value_in);
    const std::string& getName();
    void calcConstantValues(Cell& cell);
    void calculateEnergy(Cell& cell, double &energy);
    void calculateFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void updateMatrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~DM_Y_Interaction(){};
private:
    Neighbors neighbors;
    std::string name,sl_r,sl_s;
    int r,s,M;
    double value,min,max;
    double value0,value1,value2,value3;
    double z_rs;
    std::complex<double> gamma_rs;
};
}
#endif
