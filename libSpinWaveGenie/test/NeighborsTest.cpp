#define BOOST_TEST_MODULE NeighborsTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <exception>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "SpinWaveGenie/Containers/Cell.h"
#include "SpinWaveGenie/Containers/UniqueThreeVectors.h"
#include "SpinWaveGenie/Genie/Neighbors.h"

using namespace SpinWaveGenie;

Cell createMnBiCell()
{
    double SA = 2.0;
    double theta = M_PI/2.0;
    
    Cell cell;
    cell.setBasisVectors(4.2827,4.2827,6.1103,90.0,90.0,120.0);
    
    Sublattice Spin0;
    std::string name = "Spin0";
    Spin0.setName(name);
    Spin0.setType("MN2");
    Spin0.setMoment(SA,theta,0.0);
    cell.addSublattice(Spin0);
    cell.addAtom(name,0.0,0.0,0.0);
    cell.addAtom(name,0.0,0.0,0.5);
    
    return cell;
}

BOOST_AUTO_TEST_CASE( FirstNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 3.0, 3.1);
    
    BOOST_CHECK(neighborlist.size() == 2);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist-3.05515,eps);
    }
}

BOOST_AUTO_TEST_CASE( SecondNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 4.2, 4.3);
    
    BOOST_CHECK(neighborlist.size() == 6);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 4.2827,eps);
    }
}

BOOST_AUTO_TEST_CASE( ThirdNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 5.2, 5.3);
    
    BOOST_CHECK(neighborlist.size() == 12);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 5.260747,eps);
    }

}

BOOST_AUTO_TEST_CASE( FourthNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 6.0, 6.2);
    
    BOOST_CHECK(neighborlist.size() == 2);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 6.1103,eps);
    }

}

BOOST_AUTO_TEST_CASE( FifthNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 7.40, 7.43);
    
    BOOST_CHECK(neighborlist.size() == 6);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 7.41785399,eps);
    }
    
}

BOOST_AUTO_TEST_CASE( SixthNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 7.44, 7.49);
    
    BOOST_CHECK(neighborlist.size() == 12);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 7.46172134,eps);
    }
    
}

BOOST_AUTO_TEST_CASE( SeventhNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 7.9, 8.1);
    
    BOOST_CHECK(neighborlist.size() == 12);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 8.02237,eps);
    }
}

BOOST_AUTO_TEST_CASE( EigthNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 8.45, 8.65);
    
    BOOST_CHECK(neighborlist.size() == 6);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 8.5654,eps);
    }

}

BOOST_AUTO_TEST_CASE( NinethNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 9.0, 9.1);
    
    BOOST_CHECK(neighborlist.size() == 12);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 9.093955,eps);
    }

}

BOOST_AUTO_TEST_CASE( TenthNeighbors )
{
    double eps = 1.0e-5;
    Cell cell = createMnBiCell();
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 9.1, 9.2);
    
    BOOST_CHECK(neighborlist.size() == 2);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        BOOST_CHECK_SMALL(dist - 9.16545,eps);
    }
}