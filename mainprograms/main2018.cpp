
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

int main() {
    //    VecDouble vec1;
    //    ReadGmsh read;
    //    GeoMesh gmesh;
    //    read.Read(gmesh,"t1.msh");
    //    
    //    VTKGeoMesh::PrintGMeshVTK(&gmesh,"teste1.vtk");
    //    
    //    gmesh.Print(std::cout);

    int nel_x = 1;
    int nel_y = 1;
    int dim = 2;
    double len = 1;
    int pOrder = 2;

    GeoMesh *gmesh = CreateGeoMesh(nel_x, nel_y, dim, len);
    gmesh->Print(cout);

    CompMesh *cmesh = CreateCompMesh(gmesh, pOrder);

    Matrix EK;
    Matrix EF;

    int nelem = cmesh->GetElementVec().size();
    for (int i = 0; i < nelem; i++) {
        CompElement *cel = cmesh->GetElement(i);
        cel->CalcStiff(EK, EF);
    }
    
    std::cout << "\nMatriz de rigidez: " << std::endl;
    EK.Print();
    
   //Analysis an(cmesh);
    //an->RunSimulation();

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



            int matid = 1;
            GeoElement *gel = new GeoElementTemplate<GeomQuad> (TopolQuad, matid, gmesh, index);
            gmesh->SetElement(index, gel);
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

CompMesh *CreateCompMesh(GeoMesh *gmesh, int pOrder) {
    CompMesh *cmesh = new CompMesh;
    gmesh->SetReference(cmesh);

    IntRule *intrule = new IntRuleQuad;
    intrule->SetOrder(pOrder);

    Matrix p(2, 2);
    p(0, 0) = 1;
    p(0, 1) = 0;
    p(1, 0) = 0;
    p(1, 1) = 1;

    Poisson *perm = new Poisson;
    perm->SetPermeability(p);
    int64_t nelem = gmesh->NumElements();
    for (int i = 0; i < nelem; i++) {
        GeoElement *gel = gmesh->Element(i);
        CompElement *cel = new CompElementTemplate<ShapeQuad> (gel->GetIndex(), cmesh, gel);

        cel->SetIntRule(intrule);
        cel->SetStatement(perm);

        std::cout << "CompElement index: " << cel->GetIndex() << std::endl;
        std::cout << "Integration order: " << cel->GetIntRule()->GetOrder() << std::endl;
    }

    return cmesh;

}
