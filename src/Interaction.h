#ifndef _Interaction_Class
#define _Interaction_Class

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "Cell.h"



class Interaction
{
public:
    //Interaction();
    //virtual void Add_Interaction(double value_in, std::string sl_r_in,std::string sl_s_in, double min_in, double max_in) = 0;
    virtual void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN) = 0;
private:
};

#endif // Included_NameModel_H 