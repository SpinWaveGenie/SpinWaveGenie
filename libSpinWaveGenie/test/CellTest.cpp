#define BOOST_TEST_MODULE CellTest
#define BOOST_TEST_MAIN
#define _USE_MATH_DEFINES
#include <cmath>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <stdexcept>
#include <string>
#include "SpinWaveGenie/Containers/Cell.h"


using namespace SpinWaveGenie;

BOOST_AUTO_TEST_CASE( CubicBasisVectors )
{
    Matrix3 ActualBasisVectors,ActualReciprocalVectors;
    ActualBasisVectors <<
    2.0,0.0,0.0,
    0.0,2.0,0.0,
    0.0,0.0,2.0;
    
    ActualReciprocalVectors <<
    M_PI,0.0,0.0,
    0.0,M_PI,0.0,
    0.0,0.0,M_PI;
    
    Cell test;
    test.setBasisVectors(2.0,2.0,2.0,90.0,90.0,90.0);
    Matrix3 BasisVectors = test.getBasisVectors();
    Matrix3 diff = BasisVectors-ActualBasisVectors;
    BOOST_CHECK_SMALL(diff.norm(),1.0e-8);
    
    Matrix3 ReciprocalVectors = test.getReciprocalVectors();
    diff = ReciprocalVectors-ActualReciprocalVectors;
    BOOST_CHECK_SMALL(diff.norm(),1.0e-8);
}

BOOST_AUTO_TEST_CASE( HexagonalBasisVectors )
{
    Matrix3 ActualBasisVectors,ActualReciprocalVectors;
    ActualBasisVectors <<
    1.0,0.0,0.0,
    -0.5,0.5*sqrt(3.0),0.0,
    0.0,0.0,1.0;
    
    ActualReciprocalVectors <<
    2.0*M_PI,2.0*M_PI/sqrt(3),0.0,
    0.0,4.0*M_PI/sqrt(3),0.0,
    0.0,0.0,2.0*M_PI;
    
    Cell test;
    test.setBasisVectors(1.0,1.0,1.0,90.0,90.0,120.0);
    Matrix3 BasisVectors = test.getBasisVectors();
    Matrix3 diff = BasisVectors-ActualBasisVectors;
    BOOST_CHECK_SMALL(diff.norm(),1.0e-8);
    
    Matrix3 ReciprocalVectors = test.getReciprocalVectors();
    diff = ReciprocalVectors-ActualReciprocalVectors;
    BOOST_CHECK_SMALL(diff.norm(),1.0e-8);
    
}

BOOST_AUTO_TEST_CASE( AddSublattice )
{
    Cell SLtest;
    Sublattice test;
    test.setMoment(2.0,M_PI/2.0,M_PI);
    test.setName("SL1");
    test.setType("Fe3");

    SLtest.addSublattice(test);
    BOOST_CHECK(SLtest.size()==1);
    BOOST_CHECK_THROW(SLtest.addSublattice(test),std::invalid_argument);
}

BOOST_AUTO_TEST_CASE( SublatticePosition )
{
    Cell SLtest;
    Sublattice test;
    test.setMoment(2.0,M_PI/2.0,M_PI);
    test.setName("SL1");
    test.setType("Fe3");
    
    Sublattice test2;
    test2.setMoment(2.0,M_PI/2.0,M_PI);
    test2.setName("SL2");
    test2.setType("Fe3");
    
    SLtest.addSublattice(test);
    SLtest.addSublattice(test2);
    BOOST_CHECK(SLtest.getPosition("SL1")==0);
    BOOST_CHECK(SLtest.getPosition("SL2")==1);
}

BOOST_AUTO_TEST_CASE( AddAtom )
{
    Cell cell;
    cell.setBasisVectors(2.0,2.0,2.0,90.0,90.0,90.0);
    Sublattice test;
    test.setMoment(2.0,M_PI/2.0,M_PI);
    test.setName("SL1");
    test.setType("Fe3");
    cell.addSublattice(test);
    cell.addAtom("SL1",0.0,0.0,0.0);
    cell.addAtom("SL1",0.5,0.5,0.5);
    
    bool foundAtom1 = false;
    bool foundAtom2 = false;
    for (Sublattice::Iterator atom=cell.getSublattice("SL1").begin(); atom!=cell.getSublattice("SL1").end(); ++atom)
    {
        Vector3 atomdistance1(atom->get<0>()-0.0,atom->get<1>()-0.0,atom->get<2>()-0.0);
        Vector3 atomdistance2(atom->get<0>()-1.0,atom->get<1>()-1.0,atom->get<2>()-1.0);
        if (atomdistance1.norm() <0.01)
            foundAtom1 = true;
        if (atomdistance2.norm() <0.01)
            foundAtom2 = true;
    }
    BOOST_CHECK(foundAtom1);
    BOOST_CHECK(foundAtom2);
}

BOOST_AUTO_TEST_CASE(Iterator)
{
    Cell cell;
    Sublattice test;
    test.setMoment(2.0,M_PI/2.0,M_PI);
    test.setName("SL1");
    test.setType("Fe3");
    
    Sublattice test2;
    test2.setMoment(2.0,M_PI/2.0,M_PI);
    test2.setName("SL2");
    test2.setType("Fe3");
    
    cell.addSublattice(test);
    cell.addSublattice(test2);
    
    bool foundSL1 = false;
    bool foundSL2 = false;
    bool foundBlank = false;
    for (Cell::Iterator sl=cell.begin(); sl!=cell.end(); ++sl)
    {
        if ( sl->getName() == "SL1")
            foundSL1 = true;
        if ( sl->getName() == "SL2" )
            foundSL2 = true;
        if ( sl->getName() == "" )
            foundBlank = true;
    }
    BOOST_CHECK(foundSL1);
    BOOST_CHECK(foundSL2);
    BOOST_CHECK(!foundBlank);
}
