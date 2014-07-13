#ifndef __Interaction_H__
#define __Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Cell/Cell.h"

//! Base class for all classes describing magnetic interactions

class Neighbors;

namespace SpinWaveGenie
{

class Interaction
{
public:
    virtual std::vector<std::string> sublattices() const = 0;
    bool operator<(const Interaction& other) const;
    bool operator==(const Interaction& other) const;
    virtual void calculateEnergy(Cell& cell, double &energy) = 0;
    virtual void calcConstantValues(Cell& cell) = 0;
    virtual void calculateFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements) = 0;
    //! virtual method for adding terms to the matrix LN
    //! \param K reciprocal lattice point currently being simulated.
    //! \param cell pointer to Cell object containing magnetic ground state information
    //! \param LN matrix used to calculate spin wave frequencies and intensities
    virtual const std::string& getName() = 0;
    virtual void updateValue(double value) = 0;
    virtual void updateMatrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN) = 0;
    virtual Interaction* do_clone() const = 0;
    Interaction() = default;
    Interaction(const Interaction&) = default;
    //Interaction(Interaction&&) = default;
    //Interaction& operator=(const Interaction&) & = default;
    //Interaction& operator=(Interaction&&) & = default;
    virtual ~Interaction(){};
private:
};

inline
Interaction* new_clone( const Interaction& o)
{
    return o.do_clone();
}
}

#endif // __Interaction_H__ 
