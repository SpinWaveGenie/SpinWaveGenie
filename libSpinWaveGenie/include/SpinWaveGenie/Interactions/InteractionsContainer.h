#ifndef __SpinWaveGenie__InteractionsContainer__
#define __SpinWaveGenie__InteractionsContainer__

#include "SpinWaveGenie/Interactions/Interaction.h"

#include <boost/iterator/indirect_iterator.hpp>

#include <iostream>
#include <vector>

namespace SpinWaveGenie
{
//! Stores results of SpinWaveGenie Calculations at the given Q-point.
/*!
 This container stores the results of a given SpinWaveGenie calculation. Each "Point"
 contains a calculated frequency and intensity. Results can be sorted by frequency,
 or filtered so that branches with identical frequencies are combined, or branches without
 significant intensity are discarded.
 */
class InteractionsContainer
{
public:
  InteractionsContainer() = default;
  InteractionsContainer(const InteractionsContainer &other);
  InteractionsContainer &operator=(const InteractionsContainer &other);
  //! Insert Point struct into container.
  void insert(std::unique_ptr<Interaction> value);
  //! \return size of InteractionsContainer container.
  std::size_t size() const { return container.size(); }
  //! Clearcontainer so that the size is zero.
  void clear();
  void sort();
  const Interaction &operator[](std::size_t bin) { return *container[bin]; }
  typedef boost::indirect_iterator<std::vector<std::unique_ptr<Interaction>>::iterator> Iterator;
  typedef boost::indirect_iterator<std::vector<std::unique_ptr<Interaction>>::const_iterator> ConstIterator;
  // typedef std::vector<std::unique_ptr<Interaction>>::iterator Iterator;
  // typedef std::vector<std::unique_ptr<Interaction>>::const_iterator ConstIterator;
  //! \return Returns an Iterator pointing to the first element in the container.
  Iterator begin() { return container.begin(); }
  //! \return Returns an Iterator pointing to the end of the container.
  Iterator end() { return container.end(); }
  //! \return Returns an ConstIterator pointing to the first element in the container.
  ConstIterator begin() const { return container.cbegin(); }
  //! \return Returns an ConstIterator pointing to the end of the container.
  ConstIterator end() const { return container.cend(); }
  //! \return Returns an ConstIterator pointing to the first element in the container.
  ConstIterator cbegin() const { return container.cbegin(); }
  //! \return Returns an ConstIterator pointing to the end of the container.
  ConstIterator cend() const { return container.cend(); }
  friend std::ostream &operator<<(std::ostream &output, const InteractionsContainer &n);

private:
  std::vector<std::unique_ptr<Interaction>> container;
};
}
#endif /* defined(__SpinWaveGEnie__InteractionsContainer__) */
