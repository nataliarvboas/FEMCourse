//
//  GeoElement.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoElement_h
#define GeoElement_h

#include "DataTypes.h"

class GMesh;

class GeoElement
{
    // geometric mesh to which the element belongs
    GeoMesh *GMesh;
    
    // material id associated with the element
    int MaterialId;
    
public:
    
    GeoElement();
    
    GeoElement(const VecInt &nodes, GeoMesh *mesh);
    
    GeoElement(const GeoElement &copy);
    
    virtual ~GeoElement();
    
    /// access methods
    
    virtual int NodeIndex(int node) = 0;
    
    virtual GeoElementSide Neighbour(int side) = 0;
    
    virtual void SetNeighbour(int side, GeoElementSide &neigh) = 0;
    
    virtual Neighbour(int side) = 0;

    /// mapping functions
    
    virtual void X(const VecDouble &xi, VecDouble &x) = 0;
    
    virtual void GradX(const VecDouble &xi, Matrix &gradx) = 0;
    
    
};
#endif /* GeoElement_h */
