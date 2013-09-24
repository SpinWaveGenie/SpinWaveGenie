#ifndef __Interaction_Class_H__
#define __Interaction_Class_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "Cell.h"


class Interaction
{
public:
    virtual void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN) = 0;
private:
};

#endif // __Interaction_Class_H__ 