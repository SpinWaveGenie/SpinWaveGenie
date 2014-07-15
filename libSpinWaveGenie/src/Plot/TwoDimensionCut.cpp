//
//  TwoDimensionCut.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 1/16/14.
//
//
#include <fstream>
#include <Eigen/Dense>
#include "SpinWaveGenie/Plot/TwoDimensionCut.h"
#include "SpinWaveGenie/Plot/EnergyResolutionFunction.h"
#include "SpinWaveGenie/Containers/Energies.h"
#include "External/ezRateProgressBar.hpp"
#include <unistd.h>
#include <boost/thread/thread.hpp>
#ifdef USE_THREADS
#include "tbb/tbb.h"
using namespace tbb;
#endif
using namespace std;

namespace SpinWaveGenie
{

using std::string; using std::vector; using std::unique_ptr;
using std::cout; using std::endl;

void TwoDimensionCut::setFilename(string name)
{
    Filename = name;
}

void TwoDimensionCut::setPlotObject(unique_ptr<SpinWavePlot> object)
{
    InstrumentResolution = move(object);
    EnergyPoints = InstrumentResolution->getEnergies().size();
}

void TwoDimensionCut::setPoints(ThreeVectors<double> pts)
{
    Kpoints.clear();
    for(auto it = pts.begin(); it!= pts.end(); it++)
    {
        Kpoints.insert(it->get<0>(),it->get<1>(),it->get<2>());
    }
}

void TwoDimensionCut::setEnergyPoints(double min, double max, size_t points)
{
    InstrumentResolution->setEnergies(Energies(min, max, points));
    EnergyPoints = points;
}


namespace Internal
{
        
#ifdef USE_THREADS
    tbb::atomic<int> counter;
#else
    int counter;
#endif
    Eigen::MatrixXd mat;
        
    void MyThread(int numberpoints)
    {
        ez::ezRateProgressBar<int> p(numberpoints);
        p.units = "Q-points";
        p.start();
        while( counter < numberpoints)
        {
            p.update(counter);
            sleep(1);
        }
        p.update(numberpoints);
    }
        
#ifdef USE_THREADS //TBB version
    class ParallelCalculations
    {
        unique_ptr<SpinWavePlot> cut;
        ThreeVectors<double> points;
        size_t energyPoints;
    public:
        ParallelCalculations(unique_ptr<SpinWavePlot> inCut, ThreeVectors<double> inPoints) : cut(move(inCut)), points(inPoints), energyPoints(cut->getEnergies().size())
        {
        }
        ParallelCalculations(const ParallelCalculations& other) : cut(move(other.cut->clone())), points(other.points), energyPoints(other.cut->getEnergies().size())
        {
        }
        void operator()( const blocked_range<size_t>& r ) const
        {
            unique_ptr<SpinWavePlot> cut2(move(cut->clone()));
            auto it = points.cbegin()+ r.begin();
            for(size_t i = r.begin(); i!=r.end(); ++i)
            {
                vector<double> data = cut2->getCut(it->get<0>(),it->get<1>(),it->get<2>());
                for (auto j=0;j!=energyPoints;j++)
                {
                    mat(j,i) = data[j];
                }
                it++;
                counter++;
            }
        }
    };
        
    Eigen::MatrixXd GenerateMatrix(unique_ptr<SpinWavePlot> cut, ThreeVectors<double> points)
    {
        mat.resize(cut->getEnergies().size(),points.size());
        parallel_for(blocked_range<size_t>(0,points.size()), ParallelCalculations(move(cut),points));
        return mat;
    }
        
#else //Serial version
    Eigen::MatrixXd GenerateMatrix(unique_ptr<SpinWavePlot> cut, ThreeVectors<double> points)
    {
        size_t energyPoints = cut->getEnergies().size();
        mat.resize(energyPoints,points.size());
        for(auto it = points.begin(); it != points.end(); it++)
        {
            double x = it->get<0>();
            double y = it->get<1>();
            double z = it->get<2>();
                
            vector<double> val = cut->getCut(x,y,z);
            size_t m = std::distance(points.begin(),it);
            for(int n=0;n<energyPoints;n++)
            {
                mat(n,m) = val[n];
            }
            counter++;
        }
        return mat;
    }
#endif
}

Eigen::MatrixXd TwoDimensionCut::getMatrix()
{
    boost::thread myThread(Internal::MyThread,Kpoints.size());
    Eigen::MatrixXd figure = Internal::GenerateMatrix(move(InstrumentResolution),Kpoints);
    myThread.join();
    return figure;
}

void TwoDimensionCut::save()
{
    Eigen::MatrixXd figure = this->getMatrix();
    std::ofstream file(Filename);
    if (file.is_open())
    {
        file << figure;
    }
    file << endl;
    file.close();
}

}