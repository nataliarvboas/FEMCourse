//
//  GeoNode.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoNode_h
#define GeoNode_h

#include "DataTypes.h"


class GeoNode
{
    VecDouble xco;
    
public:
    
    GeoNode(VecDouble &co) : xco(co)
    {
        
    }
    
    GeoNode(const GeoNode &copy) : xco(copy.xco)
    {
        
    }
    
    GeoNode &operator=(const GeoNode &copy)
    {
        xco = copy.xco;
        return *this;
    }
    
    VecDouble Co()
    {
        return xco;
    }
    
    double Coord(int co)
    {
        return xco[co];
    }
    
    void SetCo(const VecDouble &co)
    {
        xco = co;
    }
    
};
#endif /* GeoNode_h */
