//
//  ShapeQuad.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef ShapeQuad_h
#define ShapeQuad_h

#include "DataTypes.h"

class ShapeQuad
{
public:
    /// computes the shape functions in function of the coordinate in parameter space and orders of the shape functions (size of orders is number of sides of the element topology)
    static void Shape(VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi);
    
    /// returns the number of shape functions associated with a side
    static int NShapeFunctions(int side, VecInt &orders);
    
    /// returns the total number of shape functions
    static int NShapeFunctions(VecInt &orders);
    
};


#endif /* ShapeQuad_h */
