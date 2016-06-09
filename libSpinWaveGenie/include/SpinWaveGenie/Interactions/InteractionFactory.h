//
//  InteractionFactory.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/27/14.
//
//

#ifndef __spin_wave_genie__InteractionFactory__
#define __spin_wave_genie__InteractionFactory__

#include "SpinWaveGenie/Containers/Matrices.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include "SpinWaveGenie/Memory.h"
#include <iostream>
#include <string>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT InteractionFactory
{
public:
  std::unique_ptr<Interaction> getExchange(const std::string &name, double value, const std::string &sl_r,
                                           const std::string &sl_s, double min, double max);
  std::unique_ptr<Interaction> getDzyaloshinskiiMoriya(const std::string &name, double value, const Vector3 &unitVector,
                                                       const std::string &sl_r, const std::string &sl_s, double min,
                                                       double max);
  std::unique_ptr<Interaction> getAnisotropy(const std::string &name, double value, const Vector3 &unitVector,
                                             const std::string &sl_r);
  std::unique_ptr<Interaction> getMagneticField(const std::string &name_in, double value_in, const Vector3 &direction,
                                                const std::string &sl_r_in);

private:
};
}
#endif /* defined(__spin_wave_genie__InteractionFactory__) */
