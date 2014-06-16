#define BOOST_TEST_MODULE NeighborsTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <exception>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "Containers/Cell.h"
#include "Containers/UniqueThreeVectors.h"
#include "Cell/Neighbors.h"

Cell createMnBiCell()
{
    double SA = 2.0;
    double theta = M_PI/2.0;
    
    Cell cell;
    cell.setBasisVectors(4.2827,4.2827,6.1103,90.0,90.0,120.0);
    
    Sublattice Spin0;
    string name = "Spin0";
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
    Neighbors neighborlist;
    neighborlist.findNeighbors(cell,"Spin0", "Spin0", 3.0, 3.1);
    
    for(auto nbr=neighborlist.begin();nbr!=neighborlist.end();nbr++)
    {
        double dist = sqrt(pow(nbr->get<0>(),2)+pow(nbr->get<1>(),2)+pow(nbr->get<2>(),2));
        cout << nbr->get<0>() << " " << nbr->get<1>() << " " << nbr->get<2>() << " " << dist << endl;
    }
    
}

BOOST_AUTO_TEST_CASE( SecondNeighbors )
{
    
    
}

BOOST_AUTO_TEST_CASE( ThirdNeighbors )
{
    
    
}

BOOST_AUTO_TEST_CASE( FourthNeighbors )
{
    
    
}

BOOST_AUTO_TEST_CASE( FifthNeighbors )
{
    
    
}

BOOST_AUTO_TEST_CASE( SixthNeighbors )
{
    
    
}