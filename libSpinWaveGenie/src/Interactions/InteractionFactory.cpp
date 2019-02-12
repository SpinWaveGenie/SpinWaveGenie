//
//  InteractionFactory.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/27/14.
//
//

#include "SpinWaveGenie/Interactions/InteractionFactory.h"
#include "SpinWaveGenie/Interactions/ExchangeInteraction.h"
#include "SpinWaveGenie/Interactions/ExchangeInteractionSameSublattice.h"
#include "SpinWaveGenie/Interactions/DM_Y_Interaction.h"
#include "SpinWaveGenie/Interactions/DM_Z_Interaction.h"
#include "SpinWaveGenie/Interactions/AnisotropyInteraction.h"
#include "SpinWaveGenie/Interactions/MagneticFieldInteraction.h"

namespace SpinWaveGenie
{

std::unique_ptr<Interaction> InteractionFactory::getExchange(const std::string &name, double value,
                                                             const std::string &sl_r, const std::string &sl_s,
                                                             double min, double max)
{
  if (sl_r == sl_s)
  {
    return std::make_unique<ExchangeInteractionSameSublattice>(name, value, sl_r, min, max);
  }
  else
  {
    return std::make_unique<ExchangeInteraction>(name, value, sl_r, sl_s, min, max);
}
}

std::unique_ptr<Interaction>
InteractionFactory::getDzyaloshinskiiMoriya(const std::string &name, double value, const Eigen::Vector3d &direction,
                                            const std::string &sl_r, const std::string &sl_s, double min, double max)
{
  double x, y, z;
  x = direction[0];
  y = direction[1];
  z = direction[2];

  if (x < 0.01 && y > 0.99 && z < 0.01)
  {
    return std::make_unique<DM_Y_Interaction>(name, value, sl_r, sl_s, min, max);
  }
  else if (x < 0.01 && y < 0.01 && z > 0.99)
  {
    return std::make_unique<DM_Z_Interaction>(name, value, sl_r, sl_s, min, max);
  }
  else
  {
    // a general DM interaction has not yet been implemented
    return std::make_unique<DM_Z_Interaction>(name, value, sl_r, sl_s, min, max);
  }
}

std::unique_ptr<Interaction> InteractionFactory::getAnisotropy(const std::string &name, double value,
                                                               const Eigen::Vector3d &direction,
                                                               const std::string &sl_r)
{
  return std::make_unique<AnisotropyInteraction>(name, value, direction, sl_r);
}

std::unique_ptr<Interaction> InteractionFactory::getAnisotropy(const std::string &name, const Eigen::Matrix3d &matrix,
                                                               const std::string &sl_r)
{
  return std::make_unique<AnisotropyInteraction>(name, matrix, sl_r);
}

std::unique_ptr<Interaction> InteractionFactory::getMagneticField(const std::string &name, double value,
                                                                  const Eigen::Vector3d &direction,
                                                                  const std::string &sl_r)
{
  return std::make_unique<MagneticFieldInteraction>(name, value, direction, sl_r);
}
}
