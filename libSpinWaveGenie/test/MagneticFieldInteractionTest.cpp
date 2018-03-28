#define BOOST_TEST_MODULE MagneticFieldInteractionTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <exception>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "SpinWaveGenie/Containers/Containers.h"
#include "SpinWaveGenie/Genie/Genie.h"
#include "SpinWaveGenie/Interactions/Interactions.h"


using namespace SpinWaveGenie;

Cell createCell()
{
    double SA = 2.0;
    double theta = M_PI/6.0; //30 degrees

    Cell cell;
    cell.setBasisVectors(1.0,1.0,1.0,90.0,90.0,90.0);

    Sublattice Spin0;
    std::string name = "Spin0";
    Spin0.setName(name);
    Spin0.setType("None");
    Spin0.setMoment(SA,theta,M_PI/2.356);
    cell.addSublattice(Spin0);
    cell.addAtom(name,0.0,0.0,0.0);

    return cell;
}

BOOST_AUTO_TEST_CASE(ClassicalEnergyTest)
{
    Cell cell(createCell());
    InteractionFactory factory;
    SpinWaveBuilder builder(cell);
    builder.addInteraction(factory.getMagneticField("H", -3.0, Eigen::Vector3d(0.0, 0.0, 1.0), "Spin0"));
    double soln = 3.0*2.0*cos(M_PI/6);
    BOOST_CHECK_CLOSE(builder.getEnergy(),soln,1.0e-5);
}

BOOST_AUTO_TEST_CASE(FirstOrderTest)
{
    Cell cell(createCell());
    InteractionFactory factory;
    SpinWaveBuilder builder(cell);
    builder.addInteraction(factory.getMagneticField("H", -3.0, Eigen::Vector3d(0.0, 0.0, 1.0), "Spin0"));
    double soln = -3.0*sin(M_PI/6.0);
    Eigen::VectorXcd result = builder.getFirstOrderTerms();
    BOOST_CHECK_CLOSE(result(0).real(),soln,1.0e-5);
    BOOST_CHECK_SMALL(result(0).imag(),1.0e-5);
    BOOST_CHECK_CLOSE(result(1).real(),soln,1.0e-5);
    BOOST_CHECK_SMALL(result(1).imag(),1.0e-5);
}

BOOST_AUTO_TEST_CASE(SecondOrderTest)
{
    Cell cell(createCell());
    InteractionFactory factory;
    std::unique_ptr<Interaction> magneticField =
        factory.getMagneticField("H", -3.0, Eigen::Vector3d(0.0, 0.0, 1.0), "Spin0");
    magneticField->calcConstantValues(cell);
    double soln = -3.0*0.5*cos(M_PI/6.0);
    Eigen::MatrixXcd result;
    result.setZero(2,2);
    magneticField->updateMatrix(Eigen::Vector3d(0.0, 0.0, 1.0), result);
    BOOST_CHECK_CLOSE(result(0,0).real(),soln,1.0e-5);
    BOOST_CHECK_SMALL(result(0,0).imag(),1.0e-5);
    BOOST_CHECK_SMALL(result(0,1).real(),1.0e-5);
    BOOST_CHECK_SMALL(result(0,1).imag(),1.0e-5);
    BOOST_CHECK_SMALL(result(1,0).real(),1.0e-5);
    BOOST_CHECK_SMALL(result(1,0).imag(),1.0e-5);
    BOOST_CHECK_CLOSE(result(1,1).real(),soln,1.0e-5);
    BOOST_CHECK_SMALL(result(1,1).imag(),1.0e-5);
}
