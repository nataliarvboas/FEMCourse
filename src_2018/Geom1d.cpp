//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//


#include "Geom1d.h"

 /// Constructor
Geom1d::Geom1d() : fNodeIndices()
{
    for (int is = 0; is<nSides; is++) {
        fNeighbours[is] = GeoElementSide();
    }
}
