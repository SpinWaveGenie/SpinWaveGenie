//
//  MagneticFieldInteraction.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 3/17/14.
//
//

#ifndef __spin_wave_genie__MagneticFieldInteraction__
#define __spin_wave_genie__MagneticFieldInteraction__

#include <iostream>
#include "Interaction.h"

class MagneticFieldInteraction: public Interaction
{
public:
    //!
    MagneticFieldInteraction(std::string name_in, double value_in,Vector3 direction, std::string sl_r_in);
    //!
    void UpdateInteraction(double value_in,Vector3 direction, std::string sl_r_in);
    //!
    void updateValue(double value_in);
    std::string getName();
    void calcConstantValues(Cell& cell);
    void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~MagneticFieldInteraction(){};
private:
    std::string name,sl_r;
    Vector3 directions;
    double value;
    int r,M;
    std::complex<double> LNrr;
};

#endif /* defined(__spin_wave_genie__MagneticFieldInteraction__) */
