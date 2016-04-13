#ifndef __SpinWaveGenie__InteractionsContainer__
#define __SpinWaveGenie__InteractionsContainer__

#include "SpinWaveGenie/Interactions/Interaction.h"

#include <boost/iterator/indirect_iterator.hpp>

#include <iostream>
#include <vector>

namespace SpinWaveGenie
{
//! Stores a collection of user-defined interactions.
/*!
 This container stores a collection of user-defined interactions
 between different sublattices. Since the interactions class is polymorphic,
 objects must be constructed on the heap and stored as pointers. This class
 is designed to hide many of these complications from other classes.
 */
class InteractionsContainer
{
public:
  InteractionsContainer() = default;
  InteractionsContainer(const InteractionsContainer &other);
  InteractionsContainer &operator=(const InteractionsContainer &other);
  InteractionsContainer(InteractionsContainer &&) = default;
  InteractionsContainer &operator=(InteractionsContainer &&) = default;
  //! \return size of the container.
  std::size_t size() const { return container.size(); }
  //! Clear container so that the size is zero.
  void clear() { container.clear(); }
  //! Insert an interaction into container.
  void insert(std::unique_ptr<Interaction> value) { container.push_back(std::move(value)); }
  //! Sort the interactions by the sublattices they interaction with.
  void sort();
  const Interaction &operator[](std::size_t pos) const { return *container[pos]; }
  // use boost indirect_iterator
  typedef boost::indirect_iterator<std::vector<std::unique_ptr<Interaction>>::iterator> Iterator;
  typedef boost::indirect_iterator<std::vector<std::unique_ptr<Interaction>>::const_iterator> ConstIterator;
  //! \return Returns an Iterator pointing to the first element in the container.
  Iterator begin() { return container.begin(); }
  //! \return Returns an Iterator pointing to the end of the container.
  Iterator end() { return container.end(); }
  //! \return Returns a ConstIterator pointing to the first element in the container.
  ConstIterator begin() const { return container.cbegin(); }
  //! \return Returns a ConstIterator pointing to the end of the container.
  ConstIterator end() const { return container.cend(); }
  //! \return Returns a ConstIterator pointing to the first element in the container.
  ConstIterator cbegin() const { return container.cbegin(); }
  //! \return Returns a ConstIterator pointing to the end of the container.
  ConstIterator cend() const { return container.cend(); }
  //! Helper function for printing the contents of an InteractionsContainer.
  friend std::ostream &operator<<(std::ostream &output, const InteractionsContainer &n);

private:
  std::vector<std::unique_ptr<Interaction>> container;
};
}
#endif /* defined(__SpinWaveGEnie__InteractionsContainer__) */
