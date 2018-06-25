//
//  GeoElementSide.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#include "GeoElement.h"

#ifndef GeoElementSide_h
#define GeoElementSide_h

class GeoElement;

class GeoElementSide
{
    // Associated element
    GeoElement *fElement;
    
    // Associated side
    int fSide;
    
public:
    
    // Default Constructor of GeoElementSide
    GeoElementSide();
    
    // Constructor of GeoElementSide
    GeoElementSide(GeoElement *element, int side) : fElement(element), fSide(side)
    {

    }
    
    // Copy constructor of GeoElementSide
    GeoElementSide(const GeoElementSide &copy);
    
    // Operator of copy 
    GeoElementSide &operator=(const GeoElementSide &copy);
    
    int operator==(const GeoElementSide &other) const {
        return fElement == other.fElement && fSide == other.fSide;
    }
    int operator!=(const GeoElementSide &other) const {
        return ! operator==(other);
    }
    
    // Return the associated element
    GeoElement *Element() const
    {
        return fElement;
    }

    // Return the associated side
    int Side() const
    {
        return fSide;
    }

    // Return neighbour element of a given side
    GeoElementSide Neighbour() const;
    
    bool DataConsistency(GeoElementSide &candidate);
    
    int Exists() const {return (fElement != 0 && fSide > -1);}
    
    // Fill in the data structure for the neighbouring information
    void SetNeighbour(const GeoElementSide &neighbour);
    
    // Verifiy if an element is a neighbour
    bool IsNeighbour(const GeoElementSide &candidate);
    
    // Define elements neighbourhood
    void IsertConnectivity(GeoElementSide &connectivity);
    
    // Vector with all Neighbours
    void AllNeighbours(std::vector<GeoElementSide> &allneigh);
    
    // Compute all corner neighbours
    void ComputeNeighbours(std::vector<GeoElementSide> &neighbour);
    
    
};
#endif /* GeoElementSide_h */