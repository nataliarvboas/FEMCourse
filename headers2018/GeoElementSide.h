//
//  GeoElementSide.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoElementSide_h
#define GeoElementSide_h

class GeoElement;

class GeoElementSide
{
    // associated element
    GeoElement *fElement;
    
    // associated side
    int fSide;
    
public:
    
    GeoElementSide();
    
    GeoElementSide(GeoElement *element, int side) : fElement(element), fSide(side)
    {
        
    }
    
    GeoElementSide(const GeoElementSide &copy);
    
    GeoElementSide &operator=(const GeoElementSide &copy);
    
    GeoElement *Element() const
    {
        return fElement;
    }
    
    int Side() const
    {
        return fSide;
    }
    
    GeoElementSide Neighbour();
    
};
#endif /* GeoElementSide_h */
