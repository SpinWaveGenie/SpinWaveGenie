#ifndef __SpinWaveBuilder_H__
#define __SpinWaveBuilder_H__ 1

#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Interactions/Interaction.h"

namespace SpinWaveGenie
{
    
class SpinWaveBuilder
{
public:
    SpinWaveBuilder();
    SpinWaveBuilder(Cell& cellIn);
    void updateCell(Cell& cellIn);
    void addInteraction(Interaction* in);
    void addInteraction(std::unique_ptr<Interaction> in);
    void updateInteraction(std::string name, double value);
    double getEnergy();
    Eigen::VectorXcd getFirstOrderTerms();
    SpinWave createElement();
private:
    Cell cell;
    boost::ptr_vector<Interaction> interactions;
};
}
#endif
