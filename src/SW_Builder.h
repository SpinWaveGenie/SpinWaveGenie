#ifndef __SW_Builder_H__
#define __SW_Builder_H__ 1

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "SpinWave.h"
#include "Interaction.h"

class SW_Builder
{
public:
    SW_Builder();
    SW_Builder(Cell& cell_in);
    void Add_Interaction(Interaction* in);
    Eigen::VectorXcd checkFirstOrderTerms();
    SpinWave Create_Element(double KX, double KY, double KZ);
private:
    Cell cell;
    boost::ptr_vector<Interaction> interactions;
    SpinWave SW;
};

#endif
