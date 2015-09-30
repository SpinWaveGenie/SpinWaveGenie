//
//  InteractionFactory.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/27/14.
//
//

#ifndef __spin_wave_genie__InteractionFactory__
#define __spin_wave_genie__InteractionFactory__

#include <iostream>
#include "SpinWaveGenie/Memory.h"
#include <string>

#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Interactions/Interaction.h"

namespace SpinWaveGenie
{

class InteractionFactory
{
public:
  std::unique_ptr<Interaction> getExchange(std::string name, double value, std::string sl_r, std::string sl_s,
                                           double min, double max);
  std::unique_ptr<Interaction> getDzyaloshinskiiMoriya(std::string name, double value, Vector3 unitVector,
                                                       std::string sl_r, std::string sl_s, double min, double max);
  std::unique_ptr<Interaction> getAnisotropy(std::string name, double value, Vector3 unitVector, std::string sl_r);
  std::unique_ptr<Interaction> getMagneticField(std::string name_in, double value_in, Vector3 direction,
                                                std::string sl_r_in);

private:
};
}
#endif /* defined(__spin_wave_genie__InteractionFactory__) */
