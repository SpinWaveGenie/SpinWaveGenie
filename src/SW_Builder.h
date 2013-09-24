#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "SpinWave.h"
#include "Interaction.h"

class SW_Builder
{
public:
    SW_Builder(boost::shared_ptr<Cell>& cell_in);
    void Add_Interaction(Interaction * in);
    SpinWave Create_Element(double KX, double KY, double KZ);
private:
    boost::shared_ptr<Cell> cell;
    boost::ptr_vector<Interaction> interactions;
    SpinWave SW;
};
