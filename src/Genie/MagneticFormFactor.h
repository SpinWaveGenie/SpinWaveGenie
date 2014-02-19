//
//  formfactor.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/20/13.
//
//

#ifndef __MagneticFormFactor__
#define __MagneticFormFactor__

#include <iostream>
#include <unordered_map>
#include <vector>

class MagneticFormFactor
{
public:
    MagneticFormFactor();
    MagneticFormFactor(std::string type);
    void setType(std::string type);
    void setType(std::vector<std::string> types, std::vector<double> weights);
    double getFormFactor(double kx, double ky, double kz);
protected:
    std::unordered_map<std::string,std::vector<double> > coefficients;
private:
    void setType(std::string type, double weight);
    void initializeMap();
    std::vector<std::vector<double> > Farray;
    std::vector<double> NormalizedWeights;
};

#endif /* defined(__MagneticFormFactor__) */
