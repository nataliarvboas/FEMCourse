//
//  VTKGeoMesh.h
//  FemSC
//
//  Created by Philippe Devloo on 17/04/18.
//

#ifndef VTKGeoMesh_h
#define VTKGeoMesh_h

#include "GeoMesh.h"
#include <string>

class VTKGeoMesh
{
    public:
    /** @brief Generate an output of all geomesh to VTK */
    static void PrintGMeshVTK(GeoMesh *gmesh, const std::string &filename);

};

#endif /* VTKGeoMesh_h */
