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

namespace SpinWaveGenie
{

//! Calculates the magnetic form factor at a given Q-point.
/*!
 This class calculates the magnetic form factor at a given Q-point. The analytical approximation \f$ F(Q) = \left< j_0 \right> \f$ is used where
 \f[
     \left< j_0 \left( s \right) \right> = Ae^{-as^2} + Be^{-bs^2} Ce^{-cs^2} + D
 \f]
 and \f$ s = \frac{Q}{4\pi} \f$. Weighted form factors containing multiple ions types can be defined for doped systems.

Types and coefficients are from page 2.5-3 of the  <a href="http://www.ill.eu/index.php?eID=tx_nawsecuredl&u=0&file=fileadmin/users_files/documents/links/documentation/NeutronDataBooklet.pdf&t=1403197249&hash=c1dbaa5d001bbf0d3fb52a2b30f81c17f34923c4"> ILL Neutron Data Booklet </a>
 */

class MagneticFormFactor
{
public:
    //! Default Constructor.
    MagneticFormFactor();
    //! Construct Magnetic Form Factor object of a given type.
    //! \param type Name of ion from the ILL Neutron Data Booklet (in CAPS)
    MagneticFormFactor(std::string type);
    MagneticFormFactor(const MagneticFormFactor& other);
    //! \Sets the type of ion.
    //! \param type Name of ion from ILL Neutron Data Booklet (in CAPS)
    void setType(std::string type);
    //! \Sets ions with multiple types.
    //! \param types Names of ions from the ILL Neutron Data Booklet (in CAPS).
    //! \param weights Weights associaed with each ion. The sum of all weights will be renormalized to 1.0.
    void setType(std::vector<std::string> types, std::vector<double> weights);
    //! Get the magnetic form factor at a specific Q-point.
    //! \param kx x-component of Q-point, in Angstroms.
    //! \param ky y-component of Q-point, in Angstroms.
    //! \param kz z-component of Q-point, in Angstroms.
    double getFormFactor(double kx, double ky, double kz);
protected:
    std::unordered_map<std::string,std::vector<double> > coefficients;
private:
    void setType(std::string type, double weight);
    void initializeMap();
    std::vector<std::vector<double> > Farray;
    std::vector<double> NormalizedWeights;
};
}
#endif /* defined(__MagneticFormFactor__) */
