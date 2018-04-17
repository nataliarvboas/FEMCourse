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
class GeoElementTemplate
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
    
    /// return the enumerated element type
    virtual ElementType Type();
};

#endif /* GeoElementTemplate_h */
