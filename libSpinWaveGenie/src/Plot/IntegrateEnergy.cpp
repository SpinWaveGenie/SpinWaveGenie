//
//  IntegrateEnergy.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 8/15/14.
//
//

#include "SpinWaveGenie/Plot/IntegrateEnergy.h"
#include "External/cubature.h"

using namespace std;

namespace SpinWaveGenie
{
    
    IntegrateEnergy::IntegrateEnergy(const IntegrateEnergy& other)
    {
        resolutionFunction = move(other.resolutionFunction->clone());
        this->centeredEnergies = other.centeredEnergies;
        this->delta = other.delta;
        this->tolerance = other.tolerance;
        this->maximumEvaluations = other.maximumEvaluations;
    }
    
    IntegrateEnergy::IntegrateEnergy(unique_ptr<SpinWavePlot> resFunction, Energies energies, double delta, double tol, int maxEvals)
    {
        this->resolutionFunction = move(resFunction);
        this->centeredEnergies = energies;
        this->delta = delta;
        this->tolerance = tol;
        this->maximumEvaluations = maxEvals;
    }
    
    int IntegrateEnergy::calculateIntegrand(unsigned dim, const double *x, unsigned fdim, double *retval)
    {
        //cout << "** " << x[0] << endl;//<< " " << x[1] << " " << x[2] << endl;
        //cout << " " << tmpx << " " << tmpy << " " << tmpz << endl;
        assert(dim==1);
        assert(fdim == std::distance(centeredEnergies.begin(),centeredEnergies.end()));
        
        Energies newEnergies;
        for(auto value = centeredEnergies.begin();value != centeredEnergies.end();value++)
        {
            newEnergies.insert(*value+*x);
        }
        
        resolutionFunction->setEnergies(newEnergies);
        
        vector<double> val = resolutionFunction->getCut(kx,ky,kz);
        
        std::copy(val.begin(),val.end(),retval);
        
        for (size_t i=0;i<val.size();i++)
        {
            if (retval[i] < 0.0)
                cout << "0.0 > fval[i] = " << retval[i] << endl;
        }
        
        return 0;
    }
    
    int IntegrateEnergy::calc(unsigned dim, const double *x, void *data, unsigned fdim, double *retval)
    {
        // Call non-static member function.
        return static_cast<IntegrateEnergy*>(data)->calculateIntegrand(dim,x,fdim,retval);
    }
    
    struct DivideValue
    {
        double value;
        DivideValue(double v)
        {
            value = 1.0/v;
        }
        void operator()(double &elem) const
        {
            elem *= value;
        }
    };
    
    
    std::vector<double> IntegrateEnergy::getCut(double kxIn, double kyIn, double kzIn)
    {
        kx = kxIn;
        ky = kyIn;
        kz = kzIn;
        
        int dim = 1;
        
        size_t energyPoints = centeredEnergies.size();
        vector<double> fval(energyPoints);
        vector<double> err(energyPoints);
        
        double minimumEnergy = -delta;
        double maximumEnergy = delta;
        
        double volume = 2.0*delta;
        
        hcubature(energyPoints,IntegrateEnergy::calc, this, dim, &minimumEnergy, &maximumEnergy, maximumEvaluations, tolerance/volume, 0, ERROR_INDIVIDUAL, &fval[0], &err[0]);
        
        //cout << "volume = " << volume << endl;
        
        for (size_t i=0;i<fval.size();i++)
        {
            if (fval[i] < 0.0)
                cout << "before: 0.0 > fval[i] = " << fval[i] << endl;
        }
        
        std::for_each(fval.begin(), fval.end(), DivideValue(volume));
        
        /*for (int i=0;i<fval.size();i++)
         {
         if (fval[i] < 0.0)
         cout << "after: 0.0 > fval[i] = " << fval[i] << endl;
         }*/
        
        return fval;
    }
    
    const Cell& IntegrateEnergy::getCell() const
    {
        return resolutionFunction->getCell();
    }
    
    const Energies& IntegrateEnergy::getEnergies()
    {
        return centeredEnergies;
    }
    
    void IntegrateEnergy::setEnergies(Energies energiesIn)
    {
        centeredEnergies = energiesIn;
    }
    
    std::unique_ptr<SpinWavePlot> IntegrateEnergy::clone()
    {
        return unique_ptr<SpinWavePlot>(new IntegrateEnergy(*this));
    }
    
}