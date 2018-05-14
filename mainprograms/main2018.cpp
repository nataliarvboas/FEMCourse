

#include <iostream>
#include <math.h>
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "ReadGmsh.h"
#include "VTKGeoMesh.h"
using std::cout;
using std::endl;
using std::cin;

int main ()
{
    VecDouble vec1;
    ReadGmsh read;
    GeoMesh gmesh;
    read.Read(gmesh,"t1.msh");
    
    VTKGeoMesh::PrintGMeshVTK(&gmesh,"teste1.vtk");
    
    gmesh.Print(std::cout);
   
    return 0;
}
