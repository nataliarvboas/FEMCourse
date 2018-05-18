//
//  IntPointData.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef IntPointData_h
#define IntPointData_h

#include "DataTypes.h"

class IntPointData
{
public:
    
    VecDouble ksi;
    
    double weight;
    
    VecDouble phi;
    
    Matrix dphidksi;
    
    VecDouble x;
    
    Matrix gradx;
    
    Matrix axes;
    
    double detjac;
    
    Matrix dphidx;
    
    VecDouble solution;
    
    Matrix dsoldksi;
    
    Matrix dsoldx;
    
};
#endif /* IntPointData_h */
