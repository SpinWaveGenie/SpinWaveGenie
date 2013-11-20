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
    double getFormFactor(double kx, double ky, double kz);
private:
    void initializeMap();
    std::unordered_map<std::string,std::vector<double> > coefficients;
    std::vector<double> F;
};

#endif /* defined(__MagneticFormFactor__) */
