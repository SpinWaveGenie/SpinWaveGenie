//
//  UniquePositions.h
//  positions
//
//  Created by Hahn, Steven E. on 1/2/14.
//
//

#ifndef __positions__UniquePositions__
#define __positions__UniquePositions__

#include <iostream>
#include "Positions.h"

class UniquePositions: public Positions
{
public:
    void insert(double x, double y, double z);
};

#endif /* defined(__positions__UniquePositions__) */
