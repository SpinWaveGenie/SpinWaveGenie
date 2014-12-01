//
//  IntegrateEnergy.cpp
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 8/15/14.
//
//

#include "SpinWaveGenie/Plot/IntegrateEnergy.h"
#include "cuba.h"

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
        setenv("CUBACORES","0",1);
        this->resolutionFunction = move(resFunction);
        this->centeredEnergies = energies;
        this->delta = delta;
        this->tolerance = tol;
        this->maximumEvaluations = maxEvals;
    }
    
    int IntegrateEnergy::calculateIntegrand(const int* ndim, const double *x, const int* ncomp, double *retval)
    {
        int dim = *ndim;
        int fdim = *ncomp;
        //cout << "** " << x[0] << endl;//<< " " << x[1] << " " << x[2] << endl;
        //assert(dim==1);
        assert(fdim == std::distance(centeredEnergies.begin(),centeredEnergies.end()));
        
        Energies newEnergies;
        for(auto value = centeredEnergies.begin();value != centeredEnergies.end();value++)
        {
            newEnergies.insert(*value+delta*(2.0*(*x)-1.0));
        }
        
        resolutionFunction->setEnergies(newEnergies);
        
        vector<double> val = resolutionFunction->getCut(kx,ky,kz);

        for(auto it = val.begin();it!=val.end();++it)
        {
            *it *= 2.0*delta;
        }
               
        std::copy(val.begin(),val.end(),retval);
        

        return 0;
    }
    
    int IntegrateEnergy::calc(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
    {
        // Call non-static member function.
        return static_cast<IntegrateEnergy*>(userdata)->calculateIntegrand(ndim,xx,ncomp,ff);
    }
    
    std::vector<double> IntegrateEnergy::getCut(double kxIn, double kyIn, double kzIn)
    {
        kx = kxIn;
        ky = kyIn;
        kz = kzIn;
        
        int dim = 1;
        
        size_t energyPoints = centeredEnergies.size();
        vector<double> fval(energyPoints);
        vector<double> error(energyPoints);
        vector<double> prob(energyPoints);
        
        int nregions;
        int neval,fail;
        Cuhre(dim,energyPoints, IntegrateEnergy::calc, this, 1, tolerance, tolerance, 0, 0, 50000,9,NULL, NULL,&nregions, &neval, &fail, &fval[0],&error[0], &prob[0]);
    
        /*Vegas(dim,energyPoints,IntegrateEnergy::calc,this,1,
              tolerance,tolerance,3, 9999,
              0 , 500000, 1000, 1000, 10000,
              1,NULL,NULL,
              &neval, &fail, &fval[0], &error[0], &prob[0]);
         */
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