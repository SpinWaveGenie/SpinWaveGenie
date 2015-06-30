#define BOOST_TEST_MODULE ResultsTest
#define BOOST_TEST_MAIN
#include <cmath>
#include <exception>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include "SpinWaveGenie/Containers/Results.h"

using namespace SpinWaveGenie;

bool areEqual(Point first, Point second)
{
    double eps = 1.0e-5;
    if (std::abs(first.frequency - second.frequency) < eps &&
        std::abs(first.intensity - second.intensity) < eps )
        return true;
    else
        return false;
}

Results getResults()
{
    Results results;
    
    for (auto i = 0;i!=10;i++)
    {
        Point pt;
        pt.frequency = (double)i;
        pt.intensity = (double)i;
        results.insert(pt);
        
        pt.frequency = (double)(i*2);
        pt.intensity = (double)(std::max(i*2-10,0));
        results.insert(pt);
    }
    return results;
}


BOOST_AUTO_TEST_CASE( iteratorAccess  )
{
    Results results = getResults();
    
    BOOST_CHECK(results.size() == 20);

    double sumFrequencies = 0.0;
    double sumIntensities = 0.0;

    for (auto it = results.begin();it!=results.end();++it)
    {
        sumFrequencies += it->frequency;
        sumIntensities += it->intensity;

    }
    
    BOOST_CHECK_CLOSE(sumFrequencies,135.0,1.0e-5);
    BOOST_CHECK_CLOSE(sumIntensities,65.0,1.0e-5);
}

BOOST_AUTO_TEST_CASE( isSorted)
{
    Results results = getResults();
    results.sort();
    
    BOOST_CHECK(results.size() == 20);
    
    double oldFrequency = -1.0; //no frequencies below zero
    for (auto it = results.begin();it!=results.end();++it)
    {
        BOOST_CHECK(oldFrequency < it->frequency);
    }
}

BOOST_AUTO_TEST_CASE( isUnique )
{
    Results results = getResults();
    BOOST_CHECK(results.size() == 20);
    results.uniqueSolutions();
    BOOST_CHECK(results.size() == 15);
    results.sort();
    BOOST_CHECK(std::adjacent_find(results.begin(),results.end(),areEqual)==results.end());
}

BOOST_AUTO_TEST_CASE( isSignificant )
{
    Results results = getResults();
    BOOST_CHECK(results.size() == 20);
    results.significantSolutions();
    BOOST_CHECK(results.size() == 13);
    
    for (auto it = results.begin();it!=results.end();++it)
    {
        BOOST_CHECK(1.0e-5 < it->frequency);
    }
}

BOOST_AUTO_TEST_CASE(PrintList)
{
    Results results = getResults();
    std::stringstream teststream;
    std::string header;
    
    teststream << results;
    std::getline(teststream,header,'\n');
    BOOST_CHECK_EQUAL("  frequency  intensity",header);
    for(auto result = results.cbegin(); result != results.cend();++result)
    {
        double frequency,intensity;
        teststream >> frequency >> intensity;
        BOOST_CHECK_CLOSE(result->frequency,frequency, 1.0e-5);
        BOOST_CHECK_CLOSE(result->intensity,intensity,1.0e-5);
    }
}

