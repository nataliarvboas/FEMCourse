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
    TGeom Geom;
    
public:
    
    GeoElementTemplate(const VecInt &nodeindices);
    
    
}
#endif /* GeoElementTemplate_h */
