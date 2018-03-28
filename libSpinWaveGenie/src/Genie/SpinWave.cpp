#include <Eigen/Cholesky>
#include <iomanip>
#include <cmath>
#include <random>

#include "SpinWaveGenie/Genie/SpinWave.h"
#include "SpinWaveGenie/Genie/Neighbors.h"

using namespace Eigen;
using namespace std;

namespace SpinWaveGenie
{

SpinWave::SpinWave(const Cell &cell_in, const InteractionsContainer &interactions_in)
    : cell(cell_in), M(cell.size()), N(2 * M), interactions(interactions_in)
{
  LN.setZero(N, N);
  SS.resize(N);
  std::fill_n(SS.data(), M, 1.0);
  std::fill_n(std::next(SS.data(), M), M, -1.0);
}

void SpinWave::createMatrix(double KX, double KY, double KZ)
{
  const Eigen::Matrix3d &recip = cell.getReciprocalVectors();
  Eigen::Vector3d K(KX, KY, KZ);
  K = K.transpose() * recip;
  this->KXP = KX;
  this->KYP = KY;
  this->KZP = KZ;
  clearMatrix();
  for (const auto &iter : interactions) // r
  {
    iter.updateMatrix(K, LN);
  }
}

void SpinWave::clearMatrix()
{
  LN.setZero();
  VI.clear();
}

const Cell &SpinWave::getCell() const { return cell; }

void SpinWave::calculateEigenvalues()
{
  // cout << LN << endl;

  LN.block(0, M, M, M) *= -1.0;
  LN.block(M, M, M, M) *= -1.0;

  LN = LN * 2.0;

  // silly fix to remove noise.
  for (size_t i = 0; i < 2 * M; i++)
  {
    for (size_t j = 0; j < 2 * M; j++)
    {
      if (std::abs(LN(i, j)) < 1.0e-10)
      {
        LN(i, j) = complex<double>(0.0, 0.0);
      }
    }
  }

  // cout << LN << endl << endl;

  ces.compute(LN);
  if (ces.info() != Success)
  {
    cout << ces.info() << endl;
  }

  for (size_t i = 0; i < N; i++)
  {
    complex<double> lambda = ces.eigenvalues()[i];
    VectorXcd v = ces.eigenvectors().col(i);
    VectorXcd tmp = lambda * v - LN * v;
    complex<double> zero_test = tmp.dot(tmp);
    if (abs(zero_test) > 1.0e-6)
    {
      cout << "Consider the first eigenvalue, lambda = " << lambda << endl;
      cout << "i = " << i << ", eigenvalue condition = " << tmp.dot(tmp) << endl;
    }
  }
}

void SpinWave::calculateWeights()
{
  using MatrixXcdRowMajor = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  MatrixXcdRowMajor XX;
  WW.resize(N);
  XX = ces.eigenvectors().adjoint();

  MatrixXi IPR(N, N);
  VectorXcd TEST(N);
  int IR;
  int IFL;
  // IPR.setZero();

  int maxIterations = 50;
  for (int ito = 0; ito < maxIterations; ito++)
  {
    // cout << "ito= " << ito << endl;
    // cout << "iteration # " << ito << endl;

    MatrixXcd ortho_test = XX * SS.asDiagonal() * XX.adjoint();
    IPR.setZero();
    IR = 1;
    for (size_t L1 = 0; L1 < N; L1++)
    {
      for (size_t L2 = 0; L2 < N; L2++)
      {
        if (L1 != L2 && abs(ortho_test(L1, L2)) > 1.0E-5)
        {
          IPR(L1, L2) = 1;
          IR = -1;
          // cout << "AR and AI matrices: " << L1 << " " << L2 << " " << ces.eigenvalues()[L1] << " " <<
          // ces.eigenvalues()[L2] << " " << ortho_test(L1,L2) << "\n";
        }
      }
    }
    // cout << ortho_test << endl;
    // cout << "IR= " << IR << endl;
    if (IR != -1)
    {
      TEST = ortho_test.diagonal();
      break;
    }
    for (size_t L1 = 0; L1 < N; L1++)
    {
      IFL = 0;
      for (size_t L2 = 0; L2 < N; L2++)
      {
        if (L2 > L1 && IPR(L1, L2) != 0 && IFL == 0)
        {
          for (size_t J = 0; J < N; J++)
          {
            XX(L1, J) -= XX(L2, J) * ortho_test(L1, L2) / ortho_test(L2, L2);
          }
          IFL = 1;
        }
      }
    }
    if (ito == maxIterations - 1)
    {
      TEST = ortho_test.diagonal();
      // cout << "Error calculating frequencies" << endl;
    }
  }

  vector<results> AL(N);
  for (std::size_t L1 = 0; L1 < N; L1++)
  {
    AL[L1].weight = TEST[L1].real();
    AL[L1].index = L1;
    XX.row(L1) /= sqrt(abs(AL[L1].weight));
  }

  //
  // Reorder the XX's by the weights
  //
  std::partition(AL.begin(), AL.end(), [](const results &i)
                 {
                   return i.weight > 0;
                 });

  // If two eigenvalues are approx. zero, std::partition
  // may fail to separate positive and negative weights.
  if (AL[M - 1].weight < 0.0 || AL[M].weight > 0.0)
  {
    std::sort(AL.begin(), AL.end(), [](const results &a, const results &b)
              {
                return a.weight > b.weight;
              });
  }

  for (size_t L1 = 0; L1 < N; L1++)
  {
    WW(L1) = ces.eigenvalues()[AL[L1].index].real(); // eigenvalue
  }

  // Swap rows to reflect ordering of eigenvalues.
  // The swap moves row L1 to a new position and the index must be
  // updated to reflect this.
  for (size_t L1 = 0; L1 < N; L1++)
  {
    auto oldPosition = std::find_if(AL.begin() + L1, AL.end(), [L1](const results &elem) -> bool
                                    {
                                      return L1 == elem.index;
                                    });
    if(oldPosition != AL.end())
    {
      XX.row(L1).swap(XX.row(AL[L1].index)); // eigenvector
      oldPosition->index = AL[L1].index;
    }
    else
    {
      throw std::runtime_error("Error while reordering eigenvectors. ");
    }
  }

  //
  // Evalue inverse of XY or XIN
  //

  XIN = XX.adjoint();
  XIN.block(0, M, M, M) *= -1.0;
  XIN.block(M, 0, M, M) *= -1.0;
}

void SpinWave::calculateIntensities()
{
  complex<double> XI(0.0, 1.0);
  double KX = KXP;
  double KY = KYP;
  double KZ = KZP;
  ArrayXXcd Intensities(M, 3);
  Intensities.setZero();

  std::size_t L2 = 0;
  for (const auto & elem : cell) // r
  {
    const Eigen::Matrix3d &V_r = elem.getInverseMatrix();
    double S_r = elem.getMoment();
    formFactor.setType(elem.getType());
    double ff = formFactor.getFormFactor(KX, KY, KZ);
    for (size_t L = 0; L < M; L++) // n
    {
      for (size_t L1 = 0; L1 < 3; L1++) // alpha
      {
        complex<double> Intensities_r =
            (V_r(L1, 0) - XI * V_r(L1, 1)) * XIN(L2, L + M) + (V_r(L1, 0) + XI * V_r(L1, 1)) * XIN(L2 + M, L + M);
        Intensities(L, L1) += sqrt(S_r) * ff * Intensities_r;
      }
    }
    ++L2;
  }

  Intensities *= Intensities.conjugate();
  Intensities *= 1.0 / (4.0 * M);

  VectorXd SXX = Intensities.col(0).real();
  VectorXd SYY = Intensities.col(1).real();
  VectorXd SZZ = Intensities.col(2).real();

  for (size_t i = 0; i < M; i++)
  {
    Point pt;
    pt.frequency = abs(WW[i + M]);
    double qSquared = pow(KX, 2) + pow(KY, 2) + pow(KZ, 2);
    if (qSquared < std::numeric_limits<double>::epsilon())
    {
      pt.intensity = std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
      pt.intensity =
          SXX(i) + SYY(i) + SZZ(i) - (pow(KX, 2) * SXX(i) + pow(KY, 2) * SYY(i) + pow(KZ, 2) * SZZ(i)) / qSquared;
    }
    VI.insert(pt);
  }
}

void SpinWave::calculate()
{
  this->calculateEigenvalues();
  this->calculateWeights();
  this->calculateIntensities();
}

}
