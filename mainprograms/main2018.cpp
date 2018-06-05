
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
#include "MathStatement.h"
#include "L2Projection.h"
#include "Poisson.h"
#include "Analysis.h"
#include "Assemble.h"

using std::cout;
using std::endl;
using std::cin;

GeoMesh *CreateGeoMesh(int nel_x, int nel_y, int dim, int len);
CompMesh *CreateCompMesh(GeoMesh *gmesh, int pOrder);
void ForceFunction(const VecDouble &co, VecDouble &result);

int main() {
    //    VecDouble vec1;
    //    ReadGmsh read;
    //    GeoMesh gmesh;
    //    read.Read(gmesh,"t1.msh");
    //    
    //    VTKGeoMesh::PrintGMeshVTK(&gmesh,"teste1.vtk");
    //    
    //    gmesh.Print(std::cout);

    int nel_x = 2;
    int nel_y = 1;
    int dim = 2;
    double len = 4;
    int pOrder = 1;

    GeoMesh *gmesh = CreateGeoMesh(nel_x, nel_y, dim, len);
    gmesh->Print(cout);

    CompMesh *cmesh = CreateCompMesh(gmesh, pOrder);

    Analysis an(cmesh);
    an.RunSimulation();

    return 0;
}

GeoMesh *CreateGeoMesh(int nel_x, int nel_y, int dim, int len) {
    int nnodes_x = nel_x + 1;
    int nnodes_y = nel_y + 1;
    int64_t nelem = nel_x*nel_y;

    GeoMesh * gmesh = new GeoMesh;

    gmesh->SetDimension(dim);
    gmesh->SetNumNodes(nnodes_x * nnodes_y);
    gmesh->SetNumElements(nel_x * nel_y);

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
            TopolQuad[1] = TopolQuad[0]+(nnodes_y);
            TopolQuad[2] = TopolQuad[0]+(nnodes_y) + 1;
            TopolQuad[3] = TopolQuad[0] + 1;

            int matid = 0;
            GeoElement *gel = new GeoElementTemplate<GeomQuad> (TopolQuad, matid, gmesh, index);
            gmesh->SetElement(index, gel);
        }
    }
    
    gmesh->BuildConnectivity();
    return gmesh;
}

CompMesh *CreateCompMesh(GeoMesh *gmesh, int pOrder) {
    CompMesh *cmesh = new CompMesh(gmesh);

    Matrix perm(2, 2);
    perm(0, 0) = 1;
    perm(0, 1) = 0;
    perm(1, 0) = 0;
    perm(1, 1) = 1;

    Poisson *p = new Poisson;
    p->SetPermeability(perm);    
    p->SetForceFunction(ForceFunction);

    cmesh->SetNumberMath(1);
    cmesh->SetMathStatement(0, p);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetNumberElement(gmesh->NumElements());

    cmesh->AutoBuild();
    return cmesh;
}

void ForceFunction(const VecDouble &co, VecDouble &result){
}