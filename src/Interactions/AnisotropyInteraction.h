//
//  Anisotropy_Interaction.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/25/13.
//
//

#ifndef __spin_wave_genie__AnisotropyInteraction__
#define __spin_wave_genie__AnisotropyInteraction__

#include <iostream>
#include "Interaction.h"

class AnisotropyInteraction: public Interaction
{
public:
    //!
    AnisotropyInteraction(double value_in,Vector3 direction, std::string sl_r_in);
    //!
    void UpdateInteraction(double value_in,Vector3 direction, std::string sl_r_in);
    //!
    void calcConstantValues(Cell& cell);
    void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN, int quadrant);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~AnisotropyInteraction(){};
private:
    std::string sl_r;
    Matrix3 directions;
    double value;
    int r,M;
    std::complex<double> LNrr,LNrrM,LNrMr,LNrMrM;
};

#endif /* defined(__spin_wave_genie__AnisotropyInteraction__) */
