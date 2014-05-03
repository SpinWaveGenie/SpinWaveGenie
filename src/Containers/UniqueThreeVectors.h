//
//  UniquePositions.h
//  positions
//
//  Created by Hahn, Steven E. on 1/2/14.
//
//

#ifndef __positions__UniqueThreeVectors__
#define __positions__UniqueThreeVectors__

#include <iostream>
#include "ThreeVectors.h"

template<class T>
class UniqueThreeVectors: public ThreeVectors<T>
{
public:
    void insert(T x, T y, T z);
};

template<class T>
void UniqueThreeVectors<T>::insert(T x, T y, T z)
{
    bool unique = true;
    double epsilon = 1.0e-6;
    for(auto it = ThreeVectors<T>::begin(); it!= ThreeVectors<T>::end(); it++)
    {
        if (std::abs(it->template get<0>() - x) < epsilon && std::abs(it->template get<1>() - y) < epsilon && std::abs(it->template get<2>() - z) < epsilon)
        {
            unique = false;
            break;
        }
    }
    if (unique)
    {
        //std::cout << "adding neighbor at "<< x << " " << y << " " << z << std::endl;
        ThreeVectors<T>::insert(x,y,z);
    }
}

#endif /* defined(__positions__UniquePositions__) */
