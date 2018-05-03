#ifndef __DMZ_Interaction_H__
#define __DMZ_Interaction_H__ 1

#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Export.h"
#include "SpinWaveGenie/Genie/Neighbors.h"
#include "SpinWaveGenie/Interactions/Interaction.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace SpinWaveGenie
{

class SPINWAVEGENIE_EXPORT DM_Z_Interaction : public Interaction
{
public:
  DM_Z_Interaction(const std::string &name_in, double value_in, const std::string &sl_r_in, const std::string &sl_s_in,
                   double min_in, double max_in);
  void updateInteraction(double value_in, const std::string &sl_r_in, const std::string &sl_s_in, double min_in,
                         double max_in);
  void updateValue(double value_in) override;
  const std::string &getName() const override;
  void calcConstantValues(const Cell &cell) override;
  void calculateEnergy(const Cell &cell, double &energy) override;
  void calculateFirstOrderTerms(const Cell &cell, Eigen::VectorXcd &elements) override;
  void updateMatrix(const Eigen::Vector3d &K, Eigen::MatrixXcd &LN) const override;
  std::array<std::string, 2> sublattices() const override;
  std::unique_ptr<Interaction> clone() const override;

private:
  Neighbors neighbors;
  std::string name, sl_r, sl_s;
  std::size_t r, s, M;
  double value, min, max;
  double tmp0, tmp1, tmp2, tmp3, tmp4;
  double z_rs;
};
}
#endif
