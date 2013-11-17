//
//  Atom.h
//  spin_wave_genie
//
//  Created by Hahn, Steven E. on 11/15/13.
//
//

#ifndef __spin_wave_genie__Atom__
#define __spin_wave_genie__Atom__

#include <iostream>

class Atom
{
public:
    Atom(double x,double y, double z);
    double dot(double x,double y, double z);
private:
    double x,y,z;
    
};

#endif /* defined(__spin_wave_genie__Atom__) */
