//
//  Results.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 6/1/14.
//
//

#ifndef __spin_wave_genie__Results__
#define __spin_wave_genie__Results__

#include <iostream>
#include <vector>

struct point
{
    double frequency;
    double intensity;
    bool operator<( const point& val ) const
    {
    	return frequency < val.frequency;
    }
};

class Results
{
public:
    void insert(point value);
    const std::size_t size() const;
    void sort();
    void clear();
    void uniqueSolutions();
    void significantSolutions();
    typedef std::vector<point>::iterator Iterator;
    typedef std::vector<point>::const_iterator ConstIterator;
    Iterator begin();
    Iterator end();
    ConstIterator cbegin();
    ConstIterator cend();
private:
    std::vector<point> results;
};

inline const size_t Results::size() const
{
    return results.size();
}

inline Results::Iterator Results::begin()
{
    return results.begin();
}

inline Results::Iterator Results::end()
{
    return results.end();
}

inline Results::ConstIterator Results::cbegin()
{
    return results.cbegin();
}

inline Results::ConstIterator Results::cend()
{
    return results.cend();
}

#endif /* defined(__spin_wave_genie__Results__) */
