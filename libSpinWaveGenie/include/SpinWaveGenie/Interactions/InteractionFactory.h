//
//  InteractionFactory.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/27/14.
//
//

#ifndef __spin_wave_genie__InteractionFactory__
#define __spin_wave_genie__InteractionFactory__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include <iostream>
#include <string>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT InteractionFactory
{
public:
  std::unique_ptr<Interaction> getExchange(const std::string &name, double value, const std::string &sl_r,
                                           const std::string &sl_s, double min, double max);
  std::unique_ptr<Interaction> getDzyaloshinskiiMoriya(const std::string &name, double value,
                                                       const Eigen::Vector3d &direction, const std::string &sl_r,
                                                       const std::string &sl_s, double min, double max);
  std::unique_ptr<Interaction> getAnisotropy(const std::string &name, double value, const Eigen::Vector3d &direction,
                                             const std::string &sl_r);
  std::unique_ptr<Interaction> getAnisotropy(const std::string &name, const Eigen::Matrix3d &matrix,
                                             const std::string &sl_r);
  std::unique_ptr<Interaction> getMagneticField(const std::string &name, double value, const Eigen::Vector3d &direction,
                                                const std::string &sl_r);

private:
};
}
#endif /* defined(__spin_wave_genie__InteractionFactory__) */
