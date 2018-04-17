//
//  GeoElement.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoElement_h
#define GeoElement_h

#include "DataTypes.h"
#include "GeoElementSide.h"
class GeoMesh;

class GeoElement
{
    
protected:
    // geometric mesh to which the element belongs
    GeoMesh *GMesh;
    
    // material id associated with the element
    int MaterialId;
    
public:
    
    GeoElement();
    
    GeoElement(int materialid, GeoMesh *mesh) : GMesh(mesh), MaterialId(materialid)
    {
        
    }
    
    GeoElement(const GeoElement &copy);
    
    virtual ~GeoElement();
    
    /// access methods
    
    virtual int NCornerNodes() = 0;
    
    virtual int NNodes() = 0;
    
    virtual int NodeIndex(int node) = 0;
    
    virtual GeoElementSide Neighbour(int side) = 0;
    
    virtual void SetNeighbour(int side, GeoElementSide &neigh) = 0;
    
    /// return the enumerated element type
    virtual ElementType Type() = 0;

    void SetMesh(GeoMesh *gmesh)
    {
        GMesh = gmesh;
    }
    
    int Material()
    {
        return MaterialId;
    }
    
    /// mapping functions
    
    virtual void X(const VecDouble &xi, VecDouble &x) = 0;
    
    virtual void GradX(const VecDouble &xi, Matrix &gradx) = 0;
    
    virtual void Print(std::ostream &out);
};
#endif /* GeoElement_h */
