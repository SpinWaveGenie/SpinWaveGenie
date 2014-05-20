#ifndef __SW_Builder_H__
#define __SW_Builder_H__ 1

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include "Genie/SpinWave.h"
#include "Interactions/Interaction.h"


class SW_Builder
{
public:
    SW_Builder();
    SW_Builder(Cell& cell_in);
    void addInteraction(Interaction* in);
    void addInteraction(std::unique_ptr<Interaction> in);
    double getEnergy();
    Eigen::VectorXcd getFirstOrderTerms();
    SpinWave Create_Element();
private:
    Cell cell;
    boost::ptr_vector<Interaction> interactions;
};

#endif
