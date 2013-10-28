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
    virtual std::vector<std::string> sublattices() const = 0;
    bool operator<(const Interaction& other) const;
    virtual void calcConstantValues(Cell& cell) = 0;
    virtual void calcChangingValues(Cell& cell, Eigen::Vector3d K) = 0;
    virtual void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements) = 0;
    //! virtual method for adding terms to the matrix LN
    //! \param K reciprocal lattice point currently being simulated.
    //! \param cell pointer to Cell object containing magnetic ground state information
    //! \param LN matrix used to calculate spin wave frequencies and intensities
    virtual void Update_Matrix(Eigen::Vector3d K, Cell& cell, Eigen::MatrixXcd &LN, int quadrant) = 0;
    virtual Interaction* do_clone() const = 0;
private:
};

inline
Interaction* new_clone( const Interaction& o)
{
    return o.do_clone();
}

#endif // __Interaction_H__ 