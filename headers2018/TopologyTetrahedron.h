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
    
    // Type definition of integration rule for tetrahedron elements
    typedef IntRuleTetrahedron LocIntRule;
    
    // Number of sides associated with tetrahedron elements elements
    const static int nSides = 15;
    
    // Number of corner nodes associated with tetrahedron elements
    const static int nCorners = 4;
    
    // Dimension of tetrahedron elements
    const static int Dimension = 3;
    
    // Number of nodes associated with a side
    static int NSideNodes(int side);
    
    // Local node index of a node associated with a side
    static int SideNodeIndex(int side, int node);
    
    // Return the enumerated element type
    static ElementType Type();
};


#endif /* TopologyTetrahedron_h */
