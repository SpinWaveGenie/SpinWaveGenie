#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "../Cell.h"

struct Anis_Z_Parameters
{
    std::string sl_r;
    double value;
};

class Anis_Z_Interaction
{
public:
    Anis_Z_Interaction();
    void Add_Interaction(double value, std::string sl_r);
    void Update_Matrix(Eigen::Vector3d K, boost::shared_ptr<Cell> cell, Eigen::MatrixXcd &LN);
private:
    std::vector<Anis_Z_Parameters> Anis_Z_array;
};

