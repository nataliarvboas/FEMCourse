//
//  TopologyTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef TopologyTriangle_h
#define TopologyTriangle_h

#include "IntRuleTriangle.h"

class TopologyTriangle
{
protected:
    
    typedef IntRuleTriangle LocIntRule;
    
    const int nSides = 7;
    
    const int nCorners = 3;
    
};


#endif /* TopologyTriangle_h */
