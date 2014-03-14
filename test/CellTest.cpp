#define BOOST_TEST_MODULE CellTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "Cell/Cell.h"
#include <iostream>

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
    Sublattice test;
    test.setMoment(2.0,M_PI/2.0,M_PI);
}

BOOST_AUTO_TEST_CASE( SublatticePosition )
{
    
}

BOOST_AUTO_TEST_CASE( AddAtom)
{
    
}

BOOST_AUTO_TEST_CASE(Iterator)
{
    
}
