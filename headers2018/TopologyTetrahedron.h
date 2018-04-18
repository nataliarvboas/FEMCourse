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
public:
    
    typedef IntRuleTetrahedron LocIntRule;
    
    const static int nSides = 15;
    
    const static int nCorners = 4;
    
protected:
    
    /// Number of nodes associated with a side
    static int NSideNodes(int side);
    
    /// local node index of a node associated with a side
    static int SideNodeIndex(int side, int node);
    
    /// return the enumerated element type
    static ElementType Type();
};


#endif /* TopologyTetrahedron_h */
