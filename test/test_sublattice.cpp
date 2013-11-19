#define BOOST_TEST_MODULE SublatticeTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "Cell.h"
#include <iostream>

BOOST_AUTO_TEST_CASE( constructors_test )
{
    Sublattice test;

    BOOST_CHECK_EQUAL(test.getName(),"");
    BOOST_CHECK_EQUAL(test.getType(),"None");
    BOOST_CHECK_CLOSE(test.getMoment(),1.0,0.001);
    BOOST_CHECK_SMALL(test.getTheta() - 0.0,0.001);
    BOOST_CHECK_SMALL(test.getPhi() - 0.0,0.001);
}

BOOST_AUTO_TEST_CASE( moment_test )
{
    Sublattice test;
    test.setMoment(2.0,M_PI/2.0,M_PI);
    BOOST_CHECK_CLOSE(test.getMoment(),2.0,1.0e-8);
    BOOST_CHECK_CLOSE(test.getTheta(),M_PI/2.0,1.0e-8);
    BOOST_CHECK_CLOSE(test.getPhi(),M_PI,1.0e-8);
}

BOOST_AUTO_TEST_CASE( rotation_matrices_test )
{
    Sublattice test;
    
    Matrix3 *RotationMatrix = test.getRotationMatrix();
    Matrix3 *InverseMatrix = test.getInverseMatrix();
    BOOST_CHECK((*RotationMatrix).isIdentity());
    BOOST_CHECK((*InverseMatrix).isIdentity());

    test.setMoment(2.0,M_PI/2.0,M_PI);
    RotationMatrix = test.getRotationMatrix();
    InverseMatrix = test.getInverseMatrix();
    Matrix3 IdentityTest = (*InverseMatrix)*(*RotationMatrix);
    BOOST_CHECK(IdentityTest.isIdentity());
}



