//
//  TopologyTetrahedron.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef TopologyTetrahedron_h
#define TopologyTetrahedron_h

#include "IntRuleTetrahedron.h"

class TopologyTetrahedron
{
protected:
    
    typedef IntRuleTetrahedron LocIntRule;
    
    const int nSides = 15;
    
    const int nCorners = 4;
    
};


#endif /* TopologyTetrahedron_h */
