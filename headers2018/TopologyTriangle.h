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
    
    /// Number of nodes associated with a side
    static int NSideNodes(int side);
    
    /// local node index of a node associated with a side
    static int SideNodeIndex(int side, int node);
    
    /// return the enumerated element type
    static ElementType Type();
};


#endif /* TopologyTriangle_h */
