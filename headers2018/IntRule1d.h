//
//  IntRule1d.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRule1d__
#define __FemSC__IntRule1d__

#include <cmath>
#include <stdio.h>
#include "IntRule.h"

class IntRule1d : public IntRule
{
    
public:
  
    // Default Constructor of integration rule 1D
    IntRule1d();
    
    // Constructor of integration rule 1D
    IntRule1d(int order);
    
    // Method to set polynomial order of the integration rule 1D
    virtual void SetOrder(int order);
    
    // Integration rule 1D method obtained from Numerical Recipes
    void gauleg(const double x1, const double x2, VecDouble &x, VecDouble &w);
    
};


#endif /* defined(__FemSC__TIntRule1d__) */