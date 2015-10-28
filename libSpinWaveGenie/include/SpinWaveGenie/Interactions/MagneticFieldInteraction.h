//
//  MagneticFieldInteraction.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/17/14.
//
//

#ifndef __spin_wave_genie__MagneticFieldInteraction__
#define __spin_wave_genie__MagneticFieldInteraction__

#include <iostream>
#include "SpinWaveGenie/Interactions/Interaction.h"

namespace SpinWaveGenie
{

class MagneticFieldInteraction : public Interaction
{
public:
  //!
  MagneticFieldInteraction(std::string name_in, double value_in, Vector3 direction, std::string sl_r_in);
  //!
  void updateInteraction(double value_in, Vector3 direction, std::string sl_r_in);
  //!
  void updateValue(double value_in) override;
  const std::string &getName() override;
  void calcConstantValues(Cell &cell) override;
  void calculateEnergy(Cell &cell, double &energy) override;
  void calculateFirstOrderTerms(Cell &cell, Eigen::VectorXcd &elements) override;
  void updateMatrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN) override;
  std::array<std::string, 2> sublattices() const override;
  virtual std::unique_ptr<Interaction> clone() const override;
  virtual ~MagneticFieldInteraction(){};

private:
  std::string name, sl_r;
  Vector3 directions;
  double value;
  std::size_t r, M;
  std::complex<double> LNrr;
};
}
#endif /* defined(__spin_wave_genie__MagneticFieldInteraction__) */
