//
//  TopologyQuad.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef TopologyQuad_h
#define TopologyQuad_h

#include "IntRuleQuad.h"

class TopologyQuad
{
protected:
    
    typedef IntRuleQuad LocIntRule;
    
    const int nSides = 9;
    
    const int nCorners = 4;
    
    /// Number of nodes associated with a side
    static int NSideNodes(int side);
    
    /// local node index of a node associated with a side
    static int SideNodeIndex(int side, int node);
    
    /// return the enumerated element type
    static ElementType Type();
};


#endif /* TopologyQuad_h */
