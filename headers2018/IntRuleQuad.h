//
//  IntRuleQuad.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __FemSC__IntRuleQuad__
#define __FemSC__IntRuleQuad__

#include <stdio.h>
#include "DataTypes.h"
#include "IntRule.h"

class IntRuleQuad : public IntRule
{
    

public:
  
  IntRuleQuad();
  
  IntRuleQuad(int order);
  
  virtual void SetOrder(int order);
   
  void gaulegQuad(const double x1, const double x2, VecDouble &x, VecDouble &w);

};


#endif /* defined(__FemSC__TIntRule1d__) */
