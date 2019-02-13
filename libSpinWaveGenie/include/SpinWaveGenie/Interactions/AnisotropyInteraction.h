//
//  Anisotropy_Interaction.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/25/13.
//
//

#ifndef __spin_wave_genie__AnisotropyInteraction__
#define __spin_wave_genie__AnisotropyInteraction__

#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include <iostream>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT AnisotropyInteraction : public Interaction
{
public:
  //!
  AnisotropyInteraction(const std::string &name_in, double value_in, const Eigen::Vector3d &direction_in,
                        const std::string &sl_r_in);
  //!
  AnisotropyInteraction(const std::string &name_in, const Eigen::Matrix3d &matrix_in, const std::string &sl_r_in);
  //!
  void updateInteraction(double value_in, const Eigen::Vector3d &direction_in, const std::string &sl_r_in);
  //!
  void updateInteraction(const Eigen::Matrix3d &matrix_in, const std::string &sl_r_in);
  //!
  void updateValue(double value_in) override;
  const std::string &getName() const override;
  void calcConstantValues(const Cell &cell) override;
  void calculateEnergy(const Cell &cell, double &energy) override;
  void calculateFirstOrderTerms(const Cell &cell, Eigen::VectorXcd &elements) override;
  void updateMatrix(const Eigen::Vector3d &K, Eigen::MatrixXcd &LN) const override;
  std::array<std::string, 2> sublattices() const override;
  std::unique_ptr<Interaction> clone() const override;

private:
  std::string name, sl_r;
  Eigen::Matrix3d directions;
  double value;
  std::size_t r, M;
  std::complex<double> LNrr, LNrrM, LNrMr;
  bool matrix_input;
};
}
#endif /* defined(__spin_wave_genie__AnisotropyInteraction__) */
