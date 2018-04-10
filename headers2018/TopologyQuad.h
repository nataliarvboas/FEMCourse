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
    
};


#endif /* TopologyQuad_h */
