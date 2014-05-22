//
//  HKLDirection.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 5/22/14.
//
//

#ifndef __spin_wave_genie__HKLDirection__
#define __spin_wave_genie__HKLDirection__

#include <iostream>
#include <vector>

struct Axis
{
    double v0,v1,v2,delta;
};

class HKLDirections
{
public:
    void addDirection(int direction, double delta);
    void addDirection(double v0, double v1, double v2, double delta);
    const size_t size() const;
    const Axis& operator[](std::size_t position);
    typedef std::vector<Axis>::iterator Iterator;
    typedef std::vector<Axis>::const_iterator ConstIterator;
    //! \return Returns an iterator pointing to the first Sublattice element
    Iterator begin();
    ConstIterator cbegin();
    //! \return Returns an iterator pointing to the final Sublattice element
    Iterator end();
    ConstIterator cend();
private:
    std::vector<Axis> integrationDirections;
};

inline const Axis& HKLDirections::operator[](std::size_t position)
{
    return integrationDirections[position];
}

inline const size_t HKLDirections::size() const
{
    return integrationDirections.size();
}

inline HKLDirections::Iterator HKLDirections::begin()
{
    return integrationDirections.begin();
}

inline HKLDirections::Iterator HKLDirections::end()
{
    return integrationDirections.end();
}

inline HKLDirections::ConstIterator HKLDirections::cbegin()
{
    return integrationDirections.cbegin();
}

inline HKLDirections::ConstIterator HKLDirections::cend()
{
    return integrationDirections.cend();
}

#endif /* defined(__spin_wave_genie__HKLDirection__) */
