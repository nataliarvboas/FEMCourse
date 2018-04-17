//
//  GeoElementTemplate.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoElementTemplate_h
#define GeoElementTemplate_h

#include "GeoElement.h"

template<class TGeom>
class GeoElementTemplate : public GeoElement
{
    
protected:
    
    TGeom Geom;
    
public:
    
    
    /// constructor
    GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh);
    
    GeoElementTemplate(const GeoElementTemplate &copy);
    
    GeoElementTemplate &operator=(const GeoElementTemplate &copy);
    
    GeoElement *Clone(GeoMesh *gmesh)
    {
        GeoElement *result = new GeoElementTemplate(*this);
        result->SetMesh(gmesh);
        return result;
    }

    virtual int NCornerNodes()
    {
        return TGeom::NNodes;
    }
    
    virtual int NNodes()
    {
        return Geom.NNodes();
    }
    
    virtual int NodeIndex(int node)
    {
        return Geom.NodeIndex(node);
    }
    
    virtual GeoElementSide Neighbour(int side)
    {
        return Geom.Neighbour(side);
    }
    
    virtual void SetNeighbour(int side, GeoElementSide &neigh)
    {
        Geom.SetNeighbour(side,neigh);
    }
    
    /// return the enumerated element type
    virtual ElementType Type();
    
    virtual void X(const VecDouble &xi, VecDouble &x);
    
    virtual void GradX(const VecDouble &xi, Matrix &gradx);
    
    virtual void Print(std::ostream &out);
};

#endif /* GeoElementTemplate_h */
