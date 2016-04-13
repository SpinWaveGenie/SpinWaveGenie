#define BOOST_TEST_MODULE ResultsTest
#define BOOST_TEST_MAIN
#include "SpinWaveGenie/Interactions/InteractionsContainer.h"
#include "SpinWaveGenie/Interactions/InteractionFactory.h"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <exception>
#include <iostream>

using namespace SpinWaveGenie;

InteractionsContainer getInteractions()
{
  InteractionFactory interactions;
  InteractionsContainer result;

  Vector3 xhat(1.0, 0.0, 0.0);
  result.insert(interactions.getExchange("J", -1.0, "sl0", "sl1", 0.4, 0.6));
  result.insert(interactions.getAnisotropy("D", 0.1, xhat, "sl1"));
  result.insert(interactions.getAnisotropy("D", 0.1, xhat, "sl0"));

  return result;
}

BOOST_AUTO_TEST_CASE(iteratorAccess)
{
  InteractionsContainer interactions = getInteractions();

  BOOST_CHECK_EQUAL(interactions.size(), 3);
}

BOOST_AUTO_TEST_CASE(isSorted)
{
  InteractionsContainer interactions = getInteractions();
  interactions.sort();

  BOOST_CHECK_EQUAL(interactions.size(), 3);

  auto first = interactions[0].sublattices();
  BOOST_CHECK_EQUAL(first[0], "sl0");
  BOOST_CHECK_EQUAL(first[1], "sl0");

  auto second = interactions[1].sublattices();
  BOOST_CHECK_EQUAL(second[0], "sl0");
  BOOST_CHECK_EQUAL(second[1], "sl1");

  auto third = interactions[2].sublattices();
  BOOST_CHECK_EQUAL(third[0], "sl1");
  BOOST_CHECK_EQUAL(third[1], "sl1");
}

BOOST_AUTO_TEST_CASE(PrintList)
{
  InteractionsContainer interactions = getInteractions();
  std::stringstream teststream;
  std::string header;
  teststream << interactions;
  std::getline(teststream, header, '\n');
  BOOST_CHECK_EQUAL("  Name  sl1   sl2", header);
  for (auto result = interactions.cbegin(); result != interactions.cend(); ++result)
  {
    std::string name, sl1, sl2;
    teststream >> name >> sl1 >> sl2;
    auto sl = result->sublattices();

    BOOST_CHECK_EQUAL(result->getName(), name);
    BOOST_CHECK_EQUAL(sl[0], sl1);
    BOOST_CHECK_EQUAL(sl[1], sl2);
  }
}
