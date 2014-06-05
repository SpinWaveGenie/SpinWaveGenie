#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <array>
#include <algorithm>
#include <memory>
#include <Eigen/Dense>
#include "nlopt.hpp"
#include "Genie/SpinWave.h"
#include "Genie/SpinWaveBuilder.h"
#include "Cell/Cell.h"
#include "Cell/Neighbors.h"
#include "Interactions/InteractionFactory.h"
#include "SpinWavePlot/EnergyResolutionFunction.h"
#include "SpinWavePlot/OneDimensionalFactory.h"
#include "SpinWavePlot/OneDimensionalShapes.h"
#include "SpinWavePlot/IntegrateThetaPhi.h"
#include "Containers/PointsAlongLine.h"
#include "SpinWavePlot/TwoDimensionCut.h"
#include <unistd.h>
#include "Containers/Energies.h"

using namespace std;

typedef Eigen::Matrix <bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;

class LoadXYZ
{
public:
    void loadFile(string filename);
    vector<double> getXAxis();
    vector<double> getYAxis();
    Eigen::MatrixXd getData();
    MatrixXb getMask();
protected:
    void readFile(string filename);
    vector<vector<double> > inputData;
    set<double> xaxis,yaxis;
    Eigen::MatrixXd data;
    MatrixXb mask;
};

void LoadXYZ::readFile(string filename)
{
    ifstream file;
    string output;
    file.open(filename.c_str());
    while(getline(file,output))
    {
        istringstream parser;
        parser.str(output);
        vector<double> tmp(3);
        parser >> tmp[0] >> tmp[1] >> tmp[2];
        inputData.push_back(tmp);
    }
    file.close();
}

map<double,int> makeKey(set<double> axis)
{
    map<double,int> key;
    
    int i = 0;
    for (set<double>::iterator it = axis.begin(); it!=axis.end();it++)
    {
        key.insert(make_pair(*it,i));
        i++;
    }
    
    return key;
}

void LoadXYZ::loadFile(string filename)
{
    this->readFile(filename);
    
    for (vector<vector<double> >::iterator it = inputData.begin();it!=inputData.end();it++)
    {
        xaxis.insert((*it)[0]);
        yaxis.insert((*it)[1]);
    }
    
    std::map<double,int> xkey,ykey;
    xkey = makeKey(xaxis);
    ykey = makeKey(yaxis);
    
    data.resize(xaxis.size(),yaxis.size());
    for (std::vector<vector<double> >::iterator it = inputData.begin();it!=inputData.end();it++)
    {
        data(xkey[(*it)[0]],ykey[(*it)[1]]) = (*it)[2];
    }
    
    mask.resize(xaxis.size(),yaxis.size());
    for( int i = 0; i<xaxis.size(); i++)
    {
        for( int j=0; j<yaxis.size();j++)
        {
            mask(i,j) = data(i,j) != 1.0e-20;
        }
    }
}

vector<double> LoadXYZ::getXAxis()
{
    vector<double> tmp(xaxis.size());
    std::copy(xaxis.begin(),xaxis.end(),tmp.begin());
    return tmp;
}

vector<double> LoadXYZ::getYAxis()
{
    vector<double> tmp(yaxis.size());
    std::copy(yaxis.begin(),yaxis.end(),tmp.begin());
    return tmp;
}

Eigen::MatrixXd LoadXYZ::getData()
{
    return data;
}

MatrixXb LoadXYZ::getMask()
{
    return mask;
}

double calculatenorm(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2,MatrixXb mask)
{
    double norm = 0;
    for (int i = 0; i!=mat1.rows(); i++)
    {
        for (int j = 0; j!=mat1.cols(); j++)
        {
            if (mask(i,j))
            {
                norm += pow(mat1(i,j) - mat2(i,j),2);
            }
        }
   
    }
    return norm;
}

LoadXYZ xyzfile;

double myfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty())
    {
        cout << "error: no gradient available" << endl;
    }
    
    double J,J1,J2;
    J = x[1];
    J2 = x[2];
    J1 = x[3];
    
    double SA = 1.5;
    
    Cell cell;
    cell.setBasisVectors(9.91510,8.84410,5.45950,90.0,107.5780,90.0);
    
    Sublattice Spin0;
    string name0 = "Spin0";
    Spin0.setName(name0);
    Spin0.setType("CR3");
    Spin0.setMoment(SA,17.568*M_PI/180.0,M_PI);
    cell.addSublattice(Spin0);
    cell.addAtom(name0,0.0,0.91226,0.25);
    cell.addAtom(name0,0.5,0.41226,0.25);
    
    Sublattice Spin1;
    string name1 = "Spin1";
    Spin1.setName(name1);
    Spin1.setType("CR3");
    Spin1.setMoment(SA,17.568*M_PI/180.0,M_PI);
    cell.addSublattice(Spin1);
    cell.addAtom(name1,0.0,0.08774,0.75);
    cell.addAtom(name1,0.5,0.58774,0.75);
    
    SpinWaveBuilder builder(cell);
    InteractionFactory interactionFactory;
    
    builder.addInteraction(interactionFactory.getExchange("J",J,name0,name1,3.1,3.2));
    //builder.addInteraction(interactionFactory.getExchange("J1",J1,name0,name0,6.6,6.7));
    //builder.addInteraction(interactionFactory.getExchange("J1",J1,name1,name1,6.6,6.7));
    //builder.addInteraction(interactionFactory.getExchange("J2",J2,name0,name1,5.6,5.7));
    
    //builder.addInteraction(interactionFactory.getAnisotropy("D", -0.01, Vector3(0.0,sin(17.568*M_PI/180.0),-cos(17.568*M_PI/180.0)),"Spin0"));
    //builder.addInteraction(interactionFactory.getAnisotropy("D", -0.01, Vector3(0.0,sin(17.568*M_PI/180.0),-cos(17.568*M_PI/180.0)),"Spin1"));
    
    vector<double> xaxis = xyzfile.getXAxis();
    vector<double> yaxis = xyzfile.getYAxis();
    
    
    Energies energies;
    for(auto it = yaxis.begin()+50;it!=yaxis.begin()+80;++it)
    {
        energies.insert(*it);
    }
    
    ThreeVectors<double> points;
    
    for(auto it = xaxis.begin()+10;it!=xaxis.begin()+70;++it)
    {
        points.insert(0.0,0.0,*it);
    }
    
    SpinWave SW = builder.Create_Element();
    OneDimensionalFactory factory;
    auto lorentz = factory.getLorentzian(0.2,0.000001);
    unique_ptr<SpinWavePlot> res(new EnergyResolutionFunction(move(lorentz), SW, energies));
    unique_ptr<SpinWavePlot> cut(new IntegrateThetaPhi(move(res),0.001));
    
    TwoDimensionCut twodimcut;
    twodimcut.setPlotObject(move(cut));
    twodimcut.setPoints(points);
    Eigen::MatrixXd result = twodimcut.getMatrix();
    
    MatrixXb mask = xyzfile.getMask();
    Eigen::MatrixXd data = xyzfile.getData();
    
    double output = calculatenorm(data.block(10,50,60,30)/x[0],result, mask.block(10,50,60,30));

    cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
    cout << output << endl;
    
    std::ofstream file("data_test.txt");
    if (file.is_open())
    {
        file << data.block(10,50,60,30);
    }
    file << endl;
    file.close();
    
    std::ofstream file2("measurement_test.txt");
    if (file2.is_open())
    {
        file2 << result*x[0];
    }
    file2 << endl;
    file2.close();
    
    return output;
}

int main()
{
    string filename = "/Users/svh/Documents/spin_wave_genie/build/NaCrGe2O6_1p8K.dat";
    xyzfile.loadFile(filename);
    
    //cout << result.rows() << " " << result.cols() << endl;
    
    //cout << data.norm() << endl;
    //Eigen::MatrixXd zeros(yaxis.size(),xaxis.size());
    //zeros.setZero();
    //cout << calculatenorm(zeros,data, mask) << endl;
    //cout << calculatenorm(result,zeros, mask) << endl;
    //cout << calculatenorm(data,result*1.0e-4, mask) << endl;
    
    std::vector<double> ub(4);
    ub[1] =  0.6;
    ub[3] =  0.0;
    ub[2] =  0.0;
    ub[0] =  1.0e-4;
    
    std::vector<double> lb(4);
    lb[1] =  0.3;
    lb[3] = -0.0;
    lb[2] = -0.0;
    lb[0] =  0.0;
    
    std::vector<double> x(4);
    
    double minf = 0.0;
    
    //vector<double> gradient;
    
    x[1] = 0.439941;
    x[3] = -0.0;
    x[2] = 0.0;
    x[0] = 2.62781e-05;
    
    nlopt::opt opt(nlopt::LN_SBPLX,4);
    opt.set_min_objective(myfunc, NULL);
    opt.set_upper_bounds(ub);
    opt.set_lower_bounds(lb);
    opt.set_ftol_abs(5.0e-5);
    opt.set_maxeval(1000);
    opt.optimize(x, minf);

}
