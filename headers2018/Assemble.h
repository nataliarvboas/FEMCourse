//
//  Assemble.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#ifndef Assemble_h
#define Assemble_h

#include "DataTypes.h"
#include "CompMesh.h"

// class geared towards the assembly of the global system of equations
class Assemble
{
    CompMesh *cmesh = 0;
    
public:
    
    Assemble();
    
    Assemble(CompMesh *mesh);
    
    Assemble(const Assemble &copy);
    
    Assemble &operator=(const Assemble &copy);
    
    /// Compute the total number of equations
    int64_t NEquations();
    
    /// Optimize the bandwidth of the global system of equations
    void OptimizeBandwidth();
    
    /// Compute the global stiffness matrix and right hand side
    Compute(Matrix &globmat, Matrix &rhs);
    
};
#endif /* Assemble_h */