#define BOOST_TEST_MODULE SublatticeTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <boost/test/unit_test.hpp>
#include "SpinWaveGenie/Containers/Sublattice.h"
#include "SpinWaveGenie/Containers/Cell.h"
#include <iostream>

using namespace SpinWaveGenie;

BOOST_AUTO_TEST_CASE( ConstructorsTest )
{
    Sublattice test;

    BOOST_CHECK_EQUAL(test.getName(),"");
    BOOST_CHECK_EQUAL(test.getType(),"None");
    BOOST_CHECK_EQUAL(test.size(), 0);
    BOOST_CHECK_SMALL(test.getMoment(),1.0e-8);
    BOOST_CHECK_SMALL(test.getTheta(),1.0e-8);
    BOOST_CHECK_SMALL(test.getPhi(),1.0e-8);
}

BOOST_AUTO_TEST_CASE( MomentTest )
{
    Sublattice test;
    test.setMoment(2.0,M_PI/2.0,M_PI);
    BOOST_CHECK_CLOSE(test.getMoment(),2.0,1.0e-8);
    BOOST_CHECK_CLOSE(test.getTheta(),M_PI/2.0,1.0e-8);
    BOOST_CHECK_CLOSE(test.getPhi(),M_PI,1.0e-8);
}

BOOST_AUTO_TEST_CASE( RotationMatricesTest )
{
    Sublattice test;
    
    Matrix3 RotationMatrix = test.getRotationMatrix();
    Matrix3 InverseMatrix = test.getInverseMatrix();
    BOOST_CHECK(RotationMatrix.isIdentity());
    BOOST_CHECK(InverseMatrix.isIdentity());
    
    RotationMatrix = test.getRotationMatrix();
    InverseMatrix = test.getInverseMatrix();
    Matrix3 IdentityTest = InverseMatrix*RotationMatrix;
    BOOST_CHECK(IdentityTest.isIdentity());
}

BOOST_AUTO_TEST_CASE( AtomsTest )
{
    Sublattice test;
    test.addAtom(0.0,0.0,0.0);
    test.addAtom(0.5,0.5,0.5);
    BOOST_CHECK_EQUAL(test.size(), 2);

    Vector3 atom0(0.0,0.0,0.0);
    bool FoundAtom0 = false;
    Vector3 atom1(0.25,0.25,0.25);
    bool FoundAtom1 = false;
    for (const auto &point : test)
    {
      Vector3 test0(atom0[0] - point[0], atom0[1] - point[1], atom0[2] - point[2]);
      Vector3 test1(atom1[0] - point[0], atom1[1] - point[1], atom1[2] - point[2]);
      if (test0.norm() < 0.01)
        FoundAtom0 = true;
      if (test1.norm() < 0.01)
        FoundAtom1 = true;
    }
    BOOST_CHECK(FoundAtom0);
    BOOST_CHECK(!FoundAtom1);
}




