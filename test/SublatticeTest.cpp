#define BOOST_TEST_MODULE SublatticeTest
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "Cell.h"
#include <iostream>

BOOST_AUTO_TEST_CASE( ConstructorsTest )
{
    Sublattice test;

    BOOST_CHECK_EQUAL(test.getName(),"");
    BOOST_CHECK_EQUAL(test.getType(),"None");
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
    
    Vector3 atom0(0.0,0.0,0.0);
    bool FoundAtom0 = false;
    Vector3 atom1(0.0,0.0,0.0);
    bool FoundAtom1 = false;
    for (Sublattice::Iterator it = test.begin();it !=test.end();it++)
    {
       if((atom0-(*it)).norm()<0.01)
           FoundAtom0 = true;
       if((atom1-(*it)).norm()<0.01)
           FoundAtom1 = true;
    }
    BOOST_CHECK(FoundAtom0);
    BOOST_CHECK(FoundAtom1);
}




