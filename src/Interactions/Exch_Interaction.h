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
    void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN);
private:
    std::string sl_r,sl_s;
    double value,min,max;
};

#endif
