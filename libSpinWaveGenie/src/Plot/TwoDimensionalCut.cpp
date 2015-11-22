//
//  TwoDimensionalCut.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//
#include <atomic>
#include <fstream>
#include <thread>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <Eigen/Dense>

#ifdef USE_THREADS
#include "tbb/tbb.h"
#endif

#include "SpinWaveGenie/Plot/SpinWavePlot.h"
#include "SpinWaveGenie/Plot/TwoDimensionalCut.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "External/ezRateProgressBar.hpp"

namespace
{
class progressBar
  {
  public:
    progressBar(std::size_t numberPoints) : m_numberPoints(numberPoints) { m_counter = 0; };
    void increment() { m_counter++; };
    void run()
    {
      ez::ezRateProgressBar<std::size_t> p(m_numberPoints);
      p.units = "Q-points";
      p.start();
      while (m_counter < m_numberPoints)
      {
        p.update(m_counter);
#ifdef _WIN32
        Sleep(1);
#else
        sleep(1);
#endif
      }
      p.update(m_numberPoints);
    }

  private:
    std::atomic<std::size_t> m_counter;
    std::size_t m_numberPoints;
  };
}

namespace SpinWaveGenie
{

void TwoDimensionalCut::setFilename(std::string name) { this->filename = name; }

void TwoDimensionalCut::setPlotObject(std::unique_ptr<SpinWavePlot> object) { this->cut = move(object); }

void TwoDimensionalCut::setPoints(ThreeVectors<double> pts) { this->points = std::move(pts); }

void TwoDimensionalCut::setEnergyPoints(double min, double max, std::size_t points)
{
  this->cut->setEnergies(Energies(min, max, points));
}

Eigen::MatrixXd TwoDimensionalCut::getMatrix()
{
  Eigen::MatrixXd mat(cut->getEnergies().size(), points.size());
  progressBar pbar(points.size());
  std::thread myThread(&progressBar::run, &pbar);
#ifdef USE_THREADS
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, points.size()),
                    [this, &pbar, &mat](const tbb::blocked_range<std::size_t> r)
                    {
                      std::unique_ptr<SpinWaveGenie::SpinWavePlot> cutclone = cut->clone();
                      for (std::size_t m = r.begin(); m < r.end(); ++m)
                      {
                        Eigen::MatrixXd::ColXpr values = mat.col(m);
                        auto it = points.begin() + m;
                        std::vector<double> val = cutclone->getCut(it->get<0>(), it->get<1>(), it->get<2>());
#ifdef _MSC_VER
                        std::copy(val.begin(), val.end(), stdext::make_checked_array_iterator(values.data(), values.size()));
#else
                        std::copy(val.begin(), val.end(), values.data());
#endif
                        pbar.increment();
                      }
                    });
#else
  for (std::size_t m = 0; m < points.size(); ++m)
  {
    Eigen::MatrixXd::ColXpr values = mat.col(m);
    auto it = points.begin() + m;
    std::vector<double> val = cut->getCut(it->get<0>(), it->get<1>(), it->get<2>());
#ifdef _MSC_VER
    std::copy(val.begin(), val.end(), stdext::make_checked_array_iterator(values.data(), values.size()));
#else
    std::copy(val.begin(), val.end(), values.data());
#endif
    pbar.increment();
  }
#endif
  myThread.join();
  return mat;
}

void TwoDimensionalCut::save()
{
  Eigen::MatrixXd figure = getMatrix();
  std::ofstream file(this->filename + ".mat");
  if (file.is_open())
  {
    file << figure << std::endl;
  }
  file.close();

  file.open(this->filename + ".x");
  if (file.is_open())
  {
    for (auto it = points.begin(); it != points.end(); ++it)
      file << it->get<0>() << "\t" << it->get<1>() << "\t" << it->get<2>() << std::endl;
  }
  file.close();
  file.open(this->filename + ".y");
  if (file.is_open())
  {
    Energies energies = cut->getEnergies();
    for (const auto & energie : energies)
      file << (energie) << std::endl;
  }
  file.close();
}
}
