//
//  Topology1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef Topology1d_h
#define Topology1d_h

#include "IntRule1d.h"

class Topology1d
{
protected:
    
    typedef IntRule1d LocIntRule;

    const int nSides = 3;
    
    const int nCorners = 2;
    
};

#endif /* Topology1d_h */
