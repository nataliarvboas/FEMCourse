//
//  IntRuleTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRuleTetrahedron__
#define __FemSC__IntRuleTetrahedron__

#include <stdio.h>
#include "TVec.h"
#include "TMatrix.h"
#include "IntRule.h"

class IntRuleTetrahedron : public IntRule
{
  
    public:
    
    // Default Constructor of integration rule for tetrahedron elements
    IntRuleTetrahedron();
    
    // Constructor of integration rule for tetrahedron elements
    IntRuleTetrahedron(int order);

    // Method to set polynomial order of the integration rule for tetrahedro elements
    virtual void SetOrder(int order);
    
};


#endif /* defined(__FemSC__TIntRule1d__) */
