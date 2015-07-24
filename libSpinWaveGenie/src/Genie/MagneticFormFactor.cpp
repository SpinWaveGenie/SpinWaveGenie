//
//  formfactor.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/20/13.
//
//
#include <exception>
#include <stdexcept>
#include <numeric>
#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include "SpinWaveGenie/Genie/MagneticFormFactor.h"

using std::exp;
using std::pow;
using std::vector;

namespace SpinWaveGenie
{

MagneticFormFactor::MagneticFormFactor() { this->initializeMap(); }

MagneticFormFactor::MagneticFormFactor(std::string type)
{
  this->initializeMap();
  this->setType(type);
}

MagneticFormFactor::MagneticFormFactor(const MagneticFormFactor &other)
    : coefficients(other.coefficients), Farray(other.Farray), NormalizedWeights(other.NormalizedWeights)
{
}

void MagneticFormFactor::initializeMap()
{
  coefficients =
#ifdef _WIN32
      std::initializer_list<std::pair<const std::string, std::vector<double>>>
#endif
      {{"NONE", vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}},
       {"AM2", vector<double>{0.474300, 21.776100, 1.580000, 5.690200, -1.077900, 4.145100, 0.021800}},
       {"AM3", vector<double>{0.423900, 19.573900, 1.457300, 5.872200, -0.905200, 3.968200, 0.023800}},
       {"AM4", vector<double>{0.373700, 17.862499, 1.352100, 6.042600, -0.751400, 3.719900, 0.025800}},
       {"AM5", vector<double>{0.295600, 17.372499, 1.452500, 6.073400, -0.775500, 3.661900, 0.027700}},
       {"AM6", vector<double>{0.230200, 16.953300, 1.486400, 6.115900, -0.745700, 3.542600, 0.029400}},
       {"AM7", vector<double>{0.360100, 12.729900, 1.964000, 5.120300, -1.356000, 3.714200, 0.031600}},
       {"CE2", vector<double>{0.295300, 17.684601, 0.292300, 6.732900, 0.431300, 5.382700, -0.019400}},
       {"CO0", vector<double>{0.413900, 16.161600, 0.601300, 4.780500, -0.151800, 0.021000, 0.134500}},
       {"CO1", vector<double>{0.099000, 33.125198, 0.364500, 15.176800, 0.547000, 5.008100, -0.010900}},
       {"CO2", vector<double>{0.433200, 14.355300, 0.585700, 4.607700, -0.038200, 0.133800, 0.017900}},
       {"CO3", vector<double>{0.390200, 12.507800, 0.632400, 4.457400, -0.150000, 0.034300, 0.127200}},
       {"CO4", vector<double>{0.351500, 10.778500, 0.677800, 4.234300, -0.038900, 0.240900, 0.009800}},
       {"CR0", vector<double>{0.113500, 45.199001, 0.348100, 19.493099, 0.547700, 7.354200, -0.009200}},
       {"CR1", vector<double>{-0.097700, 0.047000, 0.454400, 26.005400, 0.557900, 7.489200, 0.083100}},
       {"CR2", vector<double>{1.202400, -0.005500, 0.415800, 20.547501, 0.603200, 6.956000, -1.221800}},
       {"CR3", vector<double>{-0.309400, 0.027400, 0.368000, 17.035500, 0.655900, 6.523600, 0.285600}},
       {"CR4", vector<double>{-0.232000, 0.043300, 0.310100, 14.951800, 0.718200, 6.172600, 0.204200}},
       {"CU0", vector<double>{0.090900, 34.983799, 0.408800, 11.443200, 0.512800, 3.824800, -0.012400}},
       {"CU1", vector<double>{0.074900, 34.965599, 0.414700, 11.764200, 0.523800, 3.849700, -0.012700}},
       {"CU2", vector<double>{0.023200, 34.968601, 0.402300, 11.564000, 0.588200, 3.842800, -0.013700}},
       {"CU3", vector<double>{0.003100, 34.907398, 0.358200, 10.913800, 0.653100, 3.827900, -0.014700}},
       {"CU4", vector<double>{-0.013200, 30.681700, 0.280100, 11.162600, 0.749000, 3.817200, -0.016500}},
       {"DY2", vector<double>{0.130800, 18.315500, 0.311800, 7.664500, 0.579500, 3.146900, -0.022600}},
       {"DY3", vector<double>{0.115700, 15.073200, 0.327000, 6.799100, 0.582100, 3.020200, -0.024900}},
       {"ER2", vector<double>{0.112200, 18.122299, 0.346200, 6.910600, 0.564900, 2.761400, -0.023500}},
       {"ER3", vector<double>{0.058600, 17.980200, 0.354000, 7.096400, 0.612600, 2.748200, -0.025100}},
       {"EU2", vector<double>{0.075500, 25.296000, 0.300100, 11.599300, 0.643800, 4.025200, -0.019600}},
       {"EU3", vector<double>{0.020400, 25.307800, 0.301000, 11.474400, 0.700500, 3.942000, -0.022000}},
       {"FE0", vector<double>{0.070600, 35.008499, 0.358900, 15.358300, 0.581900, 5.560600, -0.011400}},
       {"FE1", vector<double>{0.125100, 34.963299, 0.362900, 15.514400, 0.522300, 5.591400, -0.010500}},
       {"FE2", vector<double>{0.026300, 34.959702, 0.366800, 15.943500, 0.618800, 5.593500, -0.011900}},
       {"FE3", vector<double>{0.397200, 13.244200, 0.629500, 4.903400, -0.031400, 0.349600, 0.004400}},
       {"FE4", vector<double>{0.378200, 11.380000, 0.655600, 4.592000, -0.034600, 0.483300, 0.000500}},
       {"GD2", vector<double>{0.063600, 25.382299, 0.303300, 11.212500, 0.652800, 3.787700, -0.019900}},
       {"GD3", vector<double>{0.018600, 25.386700, 0.289500, 11.142100, 0.713500, 3.752000, -0.021700}},
       {"HO2", vector<double>{0.099500, 18.176100, 0.330500, 7.855600, 0.592100, 2.979900, -0.023000}},
       {"HO3", vector<double>{0.056600, 18.317600, 0.336500, 7.688000, 0.631700, 2.942700, -0.024800}},
       {"MN0", vector<double>{0.243800, 24.962900, 0.147200, 15.672800, 0.618900, 6.540300, -0.010500}},
       {"MN1", vector<double>{-0.013800, 0.421300, 0.423100, 24.667999, 0.590500, 6.654500, -0.001000}},
       {"MN2", vector<double>{0.422000, 17.684000, 0.594800, 6.005000, 0.004300, -0.609000, -0.021900}},
       {"MN3", vector<double>{0.419800, 14.282900, 0.605400, 5.468900, 0.924100, -0.008800, -0.949800}},
       {"MN4", vector<double>{0.376000, 12.566100, 0.660200, 5.132900, -0.037200, 0.563000, 0.001100}},
       {"MO0", vector<double>{0.180600, 49.056801, 1.230600, 14.785900, -0.426800, 6.986600, 0.017100}},
       {"MO1", vector<double>{0.350000, 48.035400, 1.030500, 15.060400, -0.392900, 7.479000, 0.013900}},
       {"NB0", vector<double>{0.394600, 49.229698, 1.319700, 14.821600, -0.726900, 9.615600, 0.012900}},
       {"NB1", vector<double>{0.457200, 49.918201, 1.027400, 15.725600, -0.496200, 9.157300, 0.011800}},
       {"ND2", vector<double>{0.164500, 25.045300, 0.252200, 11.978200, 0.601200, 4.946100, -0.018000}},
       {"ND3", vector<double>{0.054000, 25.029301, 0.310100, 12.102000, 0.657500, 4.722300, -0.021600}},
       {"NI0", vector<double>{-0.017200, 35.739201, 0.317400, 14.268900, 0.713600, 4.566100, -0.014300}},
       {"NI1", vector<double>{0.070500, 35.856098, 0.398400, 13.804200, 0.542700, 4.396500, -0.011800}},
       {"NI2", vector<double>{0.016300, 35.882599, 0.391600, 13.223300, 0.605200, 4.338800, -0.013300}},
       {"NI3", vector<double>{-0.013400, 35.867699, 0.267800, 12.332600, 0.761400, 4.236900, -0.016200}},
       {"NI4", vector<double>{-0.009000, 35.861401, 0.277600, 11.790400, 0.747400, 4.201100, -0.016300}},
       {"NP3", vector<double>{0.515700, 20.865400, 2.278400, 5.893000, -1.816300, 4.845700, 0.021100}},
       {"NP4", vector<double>{0.420600, 19.804600, 2.800400, 5.978300, -2.243600, 4.984800, 0.022800}},
       {"NP5", vector<double>{0.369200, 18.190001, 3.151000, 5.850000, -2.544600, 4.916400, 0.024800}},
       {"NP6", vector<double>{0.292900, 17.561100, 3.486600, 5.784700, -2.806600, 4.870700, 0.026700}},
       {"OS5", vector<double>{0.305500, 22.152000, 1.039500, 12.52900,  0.915800, 3.016000, 0.575000}},
       {"PD0", vector<double>{0.200300, 29.363300, 1.144600, 9.599300, -0.368900, 4.042300, 0.025100}},
       {"PD1", vector<double>{0.503300, 24.503700, 1.998200, 6.908200, -1.524000, 5.513300, 0.021300}},
       {"PU3", vector<double>{0.384000, 16.679300, 3.104900, 5.421000, -2.514800, 4.551200, 0.026300}},
       {"PU4", vector<double>{0.493400, 16.835501, 1.639400, 5.638400, -1.158100, 4.139900, 0.024800}},
       {"PU5", vector<double>{0.388800, 16.559200, 2.036200, 5.656700, -1.451500, 4.255200, 0.026700}},
       {"PU6", vector<double>{0.317200, 16.050699, 3.465400, 5.350700, -2.810200, 4.513300, 0.028100}},
       {"RH0", vector<double>{0.097600, 49.882500, 1.160100, 11.830700, -0.278900, 4.126600, 0.023400}},
       {"RH1", vector<double>{0.334200, 29.756399, 1.220900, 9.438400, -0.575500, 5.332000, 0.021000}},
       {"RU0", vector<double>{0.106900, 49.423801, 1.191200, 12.741700, -0.317600, 4.912500, 0.021300}},
       {"RU1", vector<double>{0.441000, 33.308601, 1.477500, 9.553100, -0.936100, 6.722000, 0.017600}},
       {"SC0", vector<double>{0.251200, 90.029602, 0.329000, 39.402100, 0.423500, 14.322200, -0.004300}},
       {"SC1", vector<double>{0.488900, 51.160301, 0.520300, 14.076400, -0.028600, 0.179200, 0.018500}},
       {"SC2", vector<double>{0.504800, 31.403500, 0.518600, 10.989700, -0.024100, 1.183100, 0.000000}},
       {"SM2", vector<double>{0.090900, 25.203199, 0.303700, 11.856200, 0.625000, 4.236600, -0.020000}},
       {"SM3", vector<double>{0.028800, 25.206800, 0.297300, 11.831100, 0.695400, 4.211700, -0.021300}},
       {"TB2", vector<double>{0.054700, 25.508600, 0.317100, 10.591100, 0.649000, 3.517100, -0.021200}},
       {"TB3", vector<double>{0.017700, 25.509501, 0.292100, 10.576900, 0.713300, 3.512200, -0.023100}},
       {"TC0", vector<double>{0.129800, 49.661098, 1.165600, 14.130700, -0.313400, 5.512900, 0.019500}},
       {"TC1", vector<double>{0.267400, 48.956600, 0.956900, 15.141300, -0.238700, 5.457800, 0.016000}},
       {"TI0", vector<double>{0.465700, 33.589802, 0.549000, 9.879100, -0.029100, 0.323200, 0.012300}},
       {"TI1", vector<double>{0.509300, 36.703300, 0.503200, 10.371300, -0.026300, 0.310600, 0.011600}},
       {"TI2", vector<double>{0.509100, 24.976299, 0.516200, 8.756900, -0.028100, 0.916000, 0.001500}},
       {"TI3", vector<double>{0.357100, 22.841299, 0.668800, 8.930600, -0.035400, 0.483300, 0.009900}},
       {"TM2", vector<double>{0.098300, 18.323601, 0.338000, 6.917800, 0.587500, 2.662200, -0.024100}},
       {"TM3", vector<double>{0.058100, 15.092200, 0.278700, 7.801500, 0.685400, 2.793100, -0.022400}},
       {"U3", vector<double>{0.505800, 23.288200, 1.346400, 7.002800, -0.872400, 4.868300, 0.019200}},
       {"U4", vector<double>{0.329100, 23.547501, 1.083600, 8.454000, -0.434000, 4.119600, 0.021400}},
       {"U5", vector<double>{0.365000, 19.803801, 3.219900, 6.281800, -2.607700, 5.301000, 0.023300}},
       {"V0", vector<double>{0.408600, 28.810900, 0.607700, 8.543700, -0.029500, 0.276800, 0.012300}},
       {"V1", vector<double>{0.444400, 32.647900, 0.568300, 9.097100, -0.228500, 0.021800, 0.215000}},
       {"V2", vector<double>{0.408500, 23.852600, 0.609100, 8.245600, -0.167600, 0.041500, 0.149600}},
       {"V3", vector<double>{0.359800, 19.336399, 0.663200, 7.617200, -0.306400, 0.029600, 0.283500}},
       {"V4", vector<double>{0.310600, 16.816000, 0.719800, 7.048700, -0.052100, 0.302000, 0.022100}},
       {"Y0", vector<double>{0.591500, 67.608101, 1.512300, 17.900400, -1.113000, 14.135900, 0.008000}},
       {"YB2", vector<double>{0.085500, 18.512300, 0.294300, 7.373400, 0.641200, 2.677700, -0.021300}},
       {"YB3", vector<double>{0.041600, 16.094900, 0.284900, 7.834100, 0.696100, 2.672500, -0.022900}},
       {"ZR0", vector<double>{0.410600, 59.996101, 1.054300, 18.647600, -0.475100, 10.540000, 0.010600}},
       {"ZR1", vector<double>{0.453200, 59.594799, 0.783400, 21.435699, -0.245100, 9.036000, 0.009800}}};
}

void MagneticFormFactor::setType(std::string type)
{
  Farray.clear();
  NormalizedWeights.clear();
  setType(type, 1.0);
}

void MagneticFormFactor::setType(std::string type, double weight)
{

  if (coefficients.find(type) != coefficients.end())
  {
    Farray.push_back(coefficients[type]);
    NormalizedWeights.push_back(weight);
  }
  else
  {
    std::cout << "Form factor coefficients for " << type << " were not found!" << std::endl;
  }
}

void MagneticFormFactor::setType(std::vector<std::string> types, std::vector<double> weights)
{

  if (types.size() != weights.size())
  {
    throw std::runtime_error("Types and weights array are not equal lengths!");
  }

  Farray.clear();
  NormalizedWeights.clear();
  double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
  if (types.size() != weights.size())
  {
    std::cout << "Size of types and weights arrays are different" << std::endl;
  }
  else
  {
    auto begin = boost::make_zip_iterator(boost::make_tuple(types.begin(), weights.begin()));
    auto end = boost::make_zip_iterator(boost::make_tuple(types.end(), weights.end()));
    for (auto element = begin; element != end; element++)
    {
      setType(element->get<0>(), element->get<1>() / sum);
    }
  }
}

class calculateFormFactor
{
public:
  calculateFormFactor(double x, double y, double z)
  {
    ms2 = -1.0 * (pow(x, 2) + pow(y, 2) + pow(z, 2)) / (16.0 * M_PI * M_PI);
  }
  double operator()(double result, const boost::tuple<vector<double>, double> &element)
  {
    const vector<double> &F = element.get<0>();
    const double weight = element.get<1>();
    return result + weight * (F[0] * exp(F[1] * ms2) + F[2] * exp(F[3] * ms2) + F[4] * exp(F[5] * ms2) + F[6]);
  }

private:
  double ms2;
};

double MagneticFormFactor::getFormFactor(double x, double y, double z)
{
  if (Farray.size() == 0)
  {
    throw std::runtime_error("Magnetic Form Factor Not Set");
  }
  auto begin = boost::make_zip_iterator(boost::make_tuple(Farray.begin(), NormalizedWeights.begin()));
  auto end = boost::make_zip_iterator(boost::make_tuple(Farray.end(), NormalizedWeights.end()));
  return std::accumulate(begin, end, 0.0, calculateFormFactor(x, y, z));
}
}
