#define BOOST_TEST_MODULE CellTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "Cell/Cell.h"
#include <iostream>
#include <stdexcept>
#include <string>


BOOST_AUTO_TEST_CASE( ConstructorsTest )
{
    
}

BOOST_AUTO_TEST_CASE( CubicBasisVectors )
{
    Matrix3 ActualBasisVectors;
    ActualBasisVectors <<
    2.0,0.0,0.0,
    0.0,2.0,0.0,
    0.0,0.0,2.0;
    
    Cell test;
    test.setBasisVectors(2.0,2.0,2.0,90.0,90.0,90.0);
    Matrix3 BasisVectors = test.getBasisVectors();
    Matrix3 diff = BasisVectors-ActualBasisVectors;
    BOOST_CHECK_SMALL(diff.norm(),1.0e-8);
}

BOOST_AUTO_TEST_CASE( HexagonalBasisVectors )
{
    Matrix3 ActualBasisVectors;
    ActualBasisVectors <<
    1.0,0.0,0.0,
    -0.5,0.5*sqrt(3.0),0.0,
    0.0,0.0,1.0;
    
    Cell test;
    test.setBasisVectors(1.0,1.0,1.0,90.0,90.0,120.0);
    Matrix3 BasisVectors = test.getBasisVectors();
    Matrix3 diff = BasisVectors-ActualBasisVectors;
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

BOOST_AUTO_TEST_CASE( AddAtom)
{
    
}

BOOST_AUTO_TEST_CASE(Iterator)
{
    
}
