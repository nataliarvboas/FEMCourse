//
//  IntRuleTriangle.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRuleTriangle__
#define __FemSC__IntRuleTriangle__

#include <stdio.h>
#include "TVec.h"
#include "DataTypes.h"
#include "IntRule.h"


class IntRuleTriangle : public IntRule
{
  
    public:
  
    // Default Constructor of integration rule for triangle elements
    IntRuleTriangle();
  
    // Constructor of integration rule for triangle elements
    IntRuleTriangle(int order);
  
    // Method to set polynomial order of the integration rule for triangle elements
    virtual void SetOrder(int order);
  
};


#endif /* defined(__FemSC__TIntRule1d__) */
