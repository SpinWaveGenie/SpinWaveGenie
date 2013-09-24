#ifndef __Interaction_H__
#define __Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "Cell.h"

//! Base class for all classes describing magnetic interactions

class Interaction
{
public:
    //! virtual method for adding terms to the matrix LN
    //! \param K reciprocal lattice point currently being simulated.
    //! \param cell pointer to Cell object containing magnetic ground state information
    //! \param LN matrix used to calculate spin wave frequencies and intensities
    virtual void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN) = 0;
private:
};

#endif // __Interaction_H__ 