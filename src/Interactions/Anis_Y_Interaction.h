#ifndef __Anis_Y_Interaction_H__
#define __Anis_Y_Interaction_H__ 1

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Cell.h"
#include "Interaction.h"

//!
/*!
 */
class Anis_Y_Interaction: public Interaction
{
public:
    //!
    Anis_Y_Interaction(double value_in, std::string sl_r_in);
    //!
    void Update_Interaction(double value_in, std::string sl_r_in);
    //!
    void calcConstantValues(Cell& cell);
    void checkFirstOrderTerms(Cell& cell, Eigen::VectorXcd &elements);
    void Update_Matrix(Eigen::Vector3d K, Eigen::MatrixXcd &LN, int quadrant);
    std::vector<std::string> sublattices() const;
    virtual Interaction* do_clone() const;
    virtual ~Anis_Y_Interaction(){};
private:
    std::string sl_r;
    double value;
    int r,M;
    std::complex<double> LNrr,LNrrM,LNrMr;
};

#endif
