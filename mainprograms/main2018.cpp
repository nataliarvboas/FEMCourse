

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
#include "GeomQuad.h"
#include "ShapeQuad.h"
#include "GeoElementTemplate.h"
#include "CompElementTemplate.h"
#include "CompElement.h"
#include "CompMesh.h"
using std::cout;
using std::endl;
using std::cin;

int main() {
    //    VecDouble vec1;
    //    ReadGmsh read;
    //    GeoMesh gmesh;
    //    read.Read(gmesh,"t1.msh");
    //    
    //    VTKGeoMesh::PrintGMeshVTK(&gmesh,"teste1.vtk");
    //    
    //    gmesh.Print(std::cout);

    int nelem_x = 1;
    int nelem_y = 1;
    int dim = 2;
    double len = 4;

    int nnodes_x = nelem_x + 1;
    int nnodes_y = nelem_y + 1;
    int64_t nelem = nelem_x*nelem_y;

    GeoMesh * gmesh = new GeoMesh;
    gmesh->SetDimension(dim);
    gmesh->SetNumNodes(nnodes_x * nnodes_y);
    gmesh->SetNumElements(nelem_x * nelem_y);

    std::vector<double> coord(3, 0.);
    int64_t id;
    for (int i = 0; i < nnodes_x; i++) {
        for (int j = 0; j < nnodes_y; j++) {
            id = i * nnodes_y + j;
            coord[0] = (i) * len / (nnodes_x - 1);
            coord[1] = (j) * len / (nnodes_y - 1);
            gmesh->Node(id).SetCo(coord);
        }
    }

    std::vector<int> TopolQuad(4, 0);
    for (int i = 0; i < (nnodes_x - 1); i++) {
        for (int j = 0; j < (nnodes_y - 1); j++) {

            int index = (i)*(nnodes_y - 1)+ (j);
            TopolQuad[0] = (i) * nnodes_y + (j);
            TopolQuad[1] = TopolQuad[0] + 1;
            TopolQuad[2] = TopolQuad[0]+(nnodes_y) + 1;
            TopolQuad[3] = TopolQuad[0]+(nnodes_y);

            int matid = 1;
            GeoElement *gel = new GeoElementTemplate<GeomQuad> (TopolQuad, matid, gmesh, index);
            gmesh->SetElement(index, gel);
        }
    }
    gmesh->BuildConnectivity();
    gmesh->Print(cout);
    

    CompMesh *cmesh = new CompMesh;
    GeoElement *gel = gmesh->Element(0);

    //GeoElement *cel = new CompElementTemplate<ShapeQuad> (0, cmesh, gel);
    

    return 0;
}
