#include <iostream>
#include <fstream>
#include <cmath>
#include <functional> 
#include <string.h>


//#include "TVec.h"
//#include "pzvec.h"
//#include "pzerror.h"
//#include "pzmanvector.h"

#include "DataTypes.h"
#include "ReadGmsh.h"
#include "VTKGeoMesh.h"

#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeoElementTemplate.h"

#include "MathStatement.h"
#include "L2Projection.h"
#include "Poisson.h"

#include "Analysis.h"
#include "Assemble.h"

#include "PostProcess.h"
#include "PostProcessTemplate.h"

using namespace std;

const double Pi = M_PI;

void ComputeConvergenceRates(std::ostream &out, VecDouble &error0, VecDouble &error, int ndiv, double l);

GeoMesh *QuadGeoMesh(int nnodes_x, int nnodes_y, double l);
GeoMesh *TriangleGeoMesh(int nnodes_x, int nnodes_y, double l);
GeoMesh *TetrahedronGeoMesh(int nnodes_x, int nnodes_y, int nnodes_z, double l);

CompMesh *CreateCompMesh(GeoMesh *gmesh, int pOrder);

void ForceFunction(const VecDouble &co, VecDouble &result);
void Sol_exact(const VecDouble &x, VecDouble &sol, Matrix &dsol);

void QuadrilateralTest(int pOrder);
void TriangleTest(int pOrder);
void TetrahedronTest(int pOrder);

int main() {
                QuadrilateralTest(1);
    //        QuadrilateralTest(2);
    //
    //        TriangleTest(1);
    //        TriangleTest(2);

    //            TetrahedronTest(1);
//    TetrahedronTest(2);


    return 0;
}

void ComputeConvergenceRates(std::ostream &out, VecDouble &error0, VecDouble &error, int ndiv, double l) {
    if (ndiv == 1) {
        error0 = error;
    }

    double h = l / ndiv;
    double h0 = 2 * h;
    out << "\n# Convergence rates #" << std::endl;
    out << "L2-Norm (u): " << (log(error[0]) - log(error0[0])) / (log(h) - log(h0)) << endl;
    out << "L2-Norm (grad u): " << (log(error[1]) - log(error0[1])) / (log(h) - log(h0)) << endl;
    out << "H1-Norm (u): " << (log(error[2]) - log(error0[2])) / (log(h) - log(h0)) << endl;

    for (int i = 0; i < error.size(); i++) {
        error0[i] = error[i];
    }
}

GeoMesh *QuadGeoMesh(int nnodes_x, int nnodes_y, double l) {
    int64_t nnodes = nnodes_x*nnodes_y;
    int64_t nelem = (nnodes_x - 1) * (nnodes_y - 1);

    GeoMesh *gmesh = new GeoMesh;

    gmesh->SetDimension(2);
    gmesh->SetNumNodes(nnodes);

    //Cria os nos da malha
    VecDouble coord(3, 0.);
    int64_t nodeid;
    for (int i = 0; i < nnodes_y; i++) {
        for (int j = 0; j < nnodes_x; j++) {
            nodeid = i * nnodes_x + j;
            coord[0] = (j) * l / (nnodes_x - 1);
            coord[1] = (i) * l / (nnodes_y - 1);
            gmesh->Node(nodeid).SetCo(coord);
        }
    }

    //Cria os elementos
    VecInt TopolQuad(4, 0);
    for (int i = 0; i < (nnodes_x - 1); i++) {
        for (int j = 0; j < (nnodes_y - 1); j++) {

            int index = i * (nnodes_y - 1) + j;
            TopolQuad[0] = i * nnodes_y + j;
            TopolQuad[1] = TopolQuad[0] + 1;
            TopolQuad[2] = TopolQuad[0]+(nnodes_y) + 1;
            TopolQuad[3] = TopolQuad[0]+(nnodes_y);

            int matid = 1;
            GeoElement *gel = new GeoElementTemplate<GeomQuad> (TopolQuad, matid, gmesh, index);
        }
    }

    VecInt topolLine(2, 0.);
    Matrix co(2, 4, 0);
    int64_t id = nelem;

    int bc0 = -1;
    int bc1 = -2;
    int bc2 = -3;
    int bc3 = -4;

    //Condicoes de contorno
    for (int64_t iel = 0; iel < nelem; iel++) {

        //Indice do no de cada corner do elemento
        int ncorners = gmesh->Element(iel)->NCornerNodes();
        for (int i = 0; i < ncorners; i++) {
            TopolQuad[i] = gmesh->Element(iel)->NodeIndex(i);
        }

        //Coordenadas x e y de cada corner do elemento
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < ncorners; j++) {
                co(i, j) = gmesh->Node(TopolQuad[j]).Co()[i];
            }
        }

        // Condicao de contorno inferior
        if (co(1, 0) == 0 && co(1, 1) == 0) {
            topolLine[0] = TopolQuad[0];
            topolLine[1] = TopolQuad[1];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc0, gmesh, id);
            id++;
        }
        // Condicao de contorno da direita
        if (co(0, 1) == l && co(0, 2) == l) {
            topolLine[0] = TopolQuad[1];
            topolLine[1] = TopolQuad[2];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc1, gmesh, id);
            id++;
        }
        // Condicao de contorno superior
        if (co(1, 2) == l && co(1, 3) == l) {
            topolLine[0] = TopolQuad[2];
            topolLine[1] = TopolQuad[3];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc2, gmesh, id);
            id++;
        }

        // Condicao de contorno da esquerda
        if (co(0, 0) == 0 && co(0, 3) == 0) {
            topolLine[0] = TopolQuad[0];
            topolLine[1] = TopolQuad[3];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc3, gmesh, id);
            id++;
        }
    }
    gmesh->BuildConnectivity();

    coord.clear();
    TopolQuad.clear();
    topolLine.clear();

    return gmesh;
}

GeoMesh *TriangleGeoMesh(int nnodes_x, int nnodes_y, double l) {
    int64_t nnodes = nnodes_x*nnodes_y;
    int64_t nelem = 2 * (nnodes_x - 1) * (nnodes_y - 1);

    GeoMesh *gmesh = new GeoMesh;

    gmesh->SetDimension(2);
    gmesh->SetNumNodes(nnodes);

    //Cria os nos da malha
    VecDouble coord(3, 0.);
    int64_t nodeid;
    for (int i = 0; i < nnodes_y; i++) {
        for (int j = 0; j < nnodes_x; j++) {
            nodeid = i * nnodes_x + j;
            coord[0] = (j) * l / (nnodes_x - 1);
            coord[1] = (i) * l / (nnodes_y - 1);
            gmesh->Node(nodeid).SetCo(coord);
        }
    }

    //Cria os elementos
    VecInt TopolTriangle_1(3, 0);
    VecInt TopolTriangle_2(3, 0);
    int index = 0;
    int matid;
    for (int i = 0; i < (nnodes_y - 1); i++) {
        for (int j = 0; j < (nnodes_x - 1); j++) {
            TopolTriangle_1[0] = i * nnodes_y + j;
            TopolTriangle_1[1] = TopolTriangle_1[0] + 1;
            TopolTriangle_1[2] = TopolTriangle_1[1] + nnodes_x;
            matid = 1;
            GeoElement *gel_1 = new GeoElementTemplate<GeomTriangle> (TopolTriangle_1, matid, gmesh, index);
            index++;

            TopolTriangle_2[0] = TopolTriangle_1[2];
            TopolTriangle_2[1] = TopolTriangle_1[0] + nnodes_x;
            TopolTriangle_2[2] = TopolTriangle_1[0];
            matid = 1;
            GeoElement *gel_2 = new GeoElementTemplate<GeomTriangle> (TopolTriangle_2, matid, gmesh, index);
            index++;
        }
    }

    VecInt topolLine(2, 0.);
    VecInt TopolTriangle(3, 0);
    Matrix co(2, 3, 0);
    int64_t id = nelem;

    int bc0 = -1;
    int bc1 = -2;
    int bc2 = -3;
    int bc3 = -4;

    //Condicoes de contorno
    for (int64_t iel = 0; iel < nelem; iel++) {

        //Indice do no de cada corner do elemento
        int ncorners = gmesh->Element(iel)->NCornerNodes();
        for (int i = 0; i < ncorners; i++) {
            TopolTriangle[i] = gmesh->Element(iel)->NodeIndex(i);
        }

        //Coordenadas x e y de cada corner do elemento
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < ncorners; j++) {
                co(i, j) = gmesh->Node(TopolTriangle[j]).Co()[i];
            }
        }

        // Condicao de contorno inferior
        if (co(1, 0) == 0 && co(1, 1) == 0) {
            topolLine[0] = TopolTriangle[0];
            topolLine[1] = TopolTriangle[1];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc0, gmesh, id);
            id++;
        }
        // Condicao de contorno da direita
        if (co(0, 1) == l && co(0, 2) == l) {
            topolLine[0] = TopolTriangle[1];
            topolLine[1] = TopolTriangle[2];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc1, gmesh, id);
            id++;
        }
        // Condicao de contorno superior
        if (co(1, 0) == l && co(1, 1) == l) {
            topolLine[0] = TopolTriangle[0];
            topolLine[1] = TopolTriangle[1];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc2, gmesh, id);
            id++;
        }

        // Condicao de contorno da esquerda
        if (co(0, 1) == 0 && co(0, 2) == 0) {
            topolLine[0] = TopolTriangle[1];
            topolLine[1] = TopolTriangle[2];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc3, gmesh, id);
            id++;
        }
    }

    gmesh->BuildConnectivity();

    coord.clear();
    TopolTriangle_1.clear();
    TopolTriangle_2.clear();
    topolLine.clear();
    TopolTriangle.clear();

    return gmesh;
}

GeoMesh *TetrahedronGeoMesh(int nnodes_x, int nnodes_y, int nnodes_z, double l) {
    int64_t nnodes = nnodes_x * nnodes_y * nnodes_z;
    int64_t nelem = (nnodes_x - 1) * (nnodes_y - 1) * (nnodes_z - 1);

    GeoMesh *gmesh = new GeoMesh;

    gmesh->SetDimension(3);
    gmesh->SetNumNodes(nnodes);

    //Cria os nos da malha
    VecDouble coord(3, 0.);
    int64_t nodeid = 0;
    for (int i = 0; i < nnodes_z; i++) {
        for (int j = 0; j < nnodes_y; j++) {
            for (int k = 0; k < nnodes_x; k++) {
                nodeid = i * nnodes_x * nnodes_y + j * nnodes_y + k;
                coord[0] = (k) * l / (nnodes_x - 1);
                coord[1] = (j) * l / (nnodes_y - 1);
                coord[2] = (i) * l / (nnodes_z - 1);
                gmesh->Node(nodeid).SetCo(coord);
            }
        }
    }

    //Cria os elementos
    VecInt TopolTetra_1(4, 0);
    VecInt TopolTetra_2(4, 0);
    VecInt TopolTetra_3(4, 0);
    VecInt TopolTetra_4(4, 0);
    VecInt TopolTetra_5(4, 0);
    VecInt TopolTetra_6(4, 0);
    int index = 0;
    int matid;
    for (int i = 0; i < (nnodes_z - 1); i++) {
        for (int j = 0; j < (nnodes_y - 1); j++) {
            for (int k = 0; k < (nnodes_x - 1); k++) {
                TopolTetra_1[0] = i * nnodes_x * nnodes_y + j * nnodes_y + k;
                TopolTetra_1[1] = TopolTetra_1[0] + 1;
                TopolTetra_1[2] = TopolTetra_1[1] + nnodes_x;
                TopolTetra_1[3] = TopolTetra_1[2] + nnodes_x*nnodes_y;
                matid = 1;
                GeoElement *gel_1 = new GeoElementTemplate<GeomTetrahedron> (TopolTetra_1, matid, gmesh, index);
                index++;

                TopolTetra_2[0] = TopolTetra_1[2];
                TopolTetra_2[1] = TopolTetra_1[0] + nnodes_x;
                TopolTetra_2[2] = TopolTetra_1[0];
                TopolTetra_2[3] = TopolTetra_1[3];
                matid = 1;
                GeoElement *gel_2 = new GeoElementTemplate<GeomTetrahedron>(TopolTetra_2, matid, gmesh, index);
                index++;

                TopolTetra_3[0] = TopolTetra_2[1];
                TopolTetra_3[1] = TopolTetra_2[3];
                TopolTetra_3[2] = TopolTetra_2[1] + nnodes_x*nnodes_y;
                TopolTetra_3[3] = TopolTetra_2[2];
                matid = 1;
                GeoElement *gel_3 = new GeoElementTemplate<GeomTetrahedron>(TopolTetra_3, matid, gmesh, index);
                index++;

                TopolTetra_4[0] = TopolTetra_3[3];
                TopolTetra_4[1] = TopolTetra_3[3] + nnodes_x*nnodes_y;
                TopolTetra_4[2] = TopolTetra_3[2];
                TopolTetra_4[3] = TopolTetra_3[1];
                matid = 1;
                GeoElement *gel_4 = new GeoElementTemplate<GeomTetrahedron>(TopolTetra_4, matid, gmesh, index);
                index++;

                TopolTetra_5[0] = TopolTetra_4[1];
                TopolTetra_5[1] = TopolTetra_4[1] + 1;
                TopolTetra_5[2] = TopolTetra_4[3];
                TopolTetra_5[3] = TopolTetra_4[0];
                matid = 1;
                GeoElement *gel_5 = new GeoElementTemplate<GeomTetrahedron>(TopolTetra_5, matid, gmesh, index);
                index++;

                TopolTetra_6[0] = TopolTetra_5[3];
                TopolTetra_6[1] = TopolTetra_5[3] + 1;
                TopolTetra_6[2] = TopolTetra_5[1];
                TopolTetra_6[3] = TopolTetra_5[2];
                matid = 1;
                GeoElement *gel_6 = new GeoElementTemplate<GeomTetrahedron>(TopolTetra_6, matid, gmesh, index);
                index++;
            }
        }
    }

    VecInt TopolTriangle(3, 0);
    VecInt TopolTetra(4, 0);
    Matrix co(3, 4, 0);
    nelem = index;
    int64_t id = nelem;

    int bc0 = -1;
    int bc1 = -2;
    int bc2 = -3;
    int bc3 = -4;
    int bc4 = -5;
    int bc5 = -6;

    //Condicoes de contorno
    for (int64_t iel = 0; iel < nelem; iel++) {

        //Indice do no de cada corner do elemento
        int ncorners = gmesh->Element(iel)->NCornerNodes();
        for (int i = 0; i < ncorners; i++) {
            TopolTetra[i] = gmesh->Element(iel)->NodeIndex(i);
        }

        //Coordenadas x e y de cada corner do elemento
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < ncorners; j++) {
                co(i, j) = gmesh->Node(TopolTetra[j]).Co()[i];
                }
        }

        //Condicao de contorno xy-inferior
        if (co(2, 0) == 0 && co(2, 1) == 0 && co(2, 2) == 0) {
            TopolTriangle[0] = TopolTetra[0];
            TopolTriangle[1] = TopolTetra[1];
            TopolTriangle[2] = TopolTetra[2];
            GeoElement *gel = new GeoElementTemplate<GeomTriangle> (TopolTriangle, bc0, gmesh, id);
            id++;
                }

        // Condicao de contorno xy-superior
        if ((co(2, 0) == l && co(2, 1) == l && co(2, 2) == l) || (co(2, 1) == l && co(2, 2) == l && co(2, 3) == l)) {
            if (co(2, 0) == l && co(2, 1) == l && co(2, 2) == l) {
                TopolTriangle[0] = TopolTetra[0];
                TopolTriangle[1] = TopolTetra[1];
                TopolTriangle[2] = TopolTetra[2];
            } else if (co(2, 1) == l && co(2, 2) == l && co(2, 3) == l) {
                TopolTriangle[0] = TopolTetra[1];
                TopolTriangle[1] = TopolTetra[2];
                TopolTriangle[2] = TopolTetra[3];
                }

            GeoElement *gel = new GeoElementTemplate<GeomTriangle> (TopolTriangle, bc1, gmesh, id);
            id++;
                }

        // Condicao de contorno yz-esquerda
        if ((co(0, 0) == 0 && co(0, 1) == 0 && co(0, 2) == 0) || (co(0, 0) == 0 && co(0, 2) == 0 && co(0, 3) == 0)) {
            if (co(0, 0) == 0 && co(0, 1) == 0 && co(0, 2) == 0) {
                TopolTriangle[0] = TopolTetra[0];
                TopolTriangle[1] = TopolTetra[1];
                TopolTriangle[2] = TopolTetra[2];
            } else if (co(0, 0) == 0 && co(0, 2) == 0 && co(0, 3) == 0) {
                TopolTriangle[0] = TopolTetra[0];
                TopolTriangle[1] = TopolTetra[2];
                TopolTriangle[2] = TopolTetra[3];
                }
            GeoElement *gel = new GeoElementTemplate<GeomTriangle> (TopolTriangle, bc2, gmesh, id);
            id++;
        }
        // Condicao de contorno yz-direita
        if (co(0, 1) == l && co(0, 2) == l && co(0, 3) == l) {
            TopolTriangle[0] = TopolTetra[1];
            TopolTriangle[1] = TopolTetra[2];
            TopolTriangle[2] = TopolTetra[3];
            GeoElement *gel = new GeoElementTemplate<GeomTriangle> (TopolTriangle, bc3, gmesh, id);
            id++;
        }

        // Condicao de contorno xz-frente
        if ((co(1, 0) == 0 && co(1, 1) == 0 && co(1, 2) == 0) || (co(1, 0) == 0 && co(1, 1) == 0 && co(1, 3) == 0)) {
            if (co(1, 0) == 0 && co(1, 1) == 0 && co(1, 2) == 0) {
                TopolTriangle[0] = TopolTetra[0];
                TopolTriangle[1] = TopolTetra[1];
                TopolTriangle[2] = TopolTetra[2];
            } else if (co(1, 0) == 0 && co(1, 1) == 0 && co(1, 3) == 0) {
                TopolTriangle[0] = TopolTetra[0];
                TopolTriangle[1] = TopolTetra[1];
                TopolTriangle[2] = TopolTetra[3];
                }
            GeoElement *gel = new GeoElementTemplate<GeomTriangle> (TopolTriangle, bc4, gmesh, id);
            id++;
            }
        // Condicao de contorno xz-tras
        if ((co(1, 0) == l && co(1, 1) == l && co(1, 2) == l) || (co(1, 0) == l && co(1, 1) == l && co(1, 3) == l)) {
            if (co(1, 0) == l && co(1, 1) == l && co(1, 2) == l) {
                TopolTriangle[0] = TopolTetra[0];
                TopolTriangle[1] = TopolTetra[1];
                TopolTriangle[2] = TopolTetra[2];
            } else if (co(1, 0) == l && co(1, 1) == l && co(1, 3) == l) {
                TopolTriangle[0] = TopolTetra[0];
                TopolTriangle[1] = TopolTetra[1];
                TopolTriangle[2] = TopolTetra[3];
        }
            GeoElement *gel = new GeoElementTemplate<GeomTriangle> (TopolTriangle, bc5, gmesh, id);
            id++;
    }
    }
    gmesh->BuildConnectivity();

    coord.clear();
    TopolTetra_1.clear();
    TopolTetra_2.clear();
    TopolTetra_3.clear();
    TopolTetra_4.clear();
    TopolTetra_5.clear();
    TopolTetra_6.clear();
    TopolTetra.clear();
    TopolTriangle.clear();

    return gmesh;

}

CompMesh *CreateCompMesh(GeoMesh *gmesh, int pOrder, int dim) {
    CompMesh *cmesh = new CompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    int nelem = cmesh->GetElementVec().size();

    Matrix perm(dim, dim);
    perm.Zero();
    for (int i = 0; i < dim; i++) {
        perm(i, i) = 1;
    }

    for (int i = 0; i < nelem; i++) {
        int matid = gmesh->Element(i)->Material();
        if (matid == 1) {
            cmesh->SetNumberMath(i + 1);
            Poisson *p = new Poisson(matid, perm);
            p->SetForceFunction(ForceFunction);
            p->SetExactSolution(Sol_exact);
            p->SetDimension(dim);
            cmesh->SetMathStatement(i, p);
        } else {
            cmesh->SetNumberMath(i + 1);
            Matrix Val1(dim, dim, 0.), Val2(dim, dim, 0.);
            L2Projection *bc = new L2Projection(0, matid, perm, Val1, Val2);
            bc->SetForceFunction(ForceFunction);
            bc->SetExactSolution(Sol_exact);
            cmesh->SetMathStatement(i, bc);
        }
    }
    cmesh->AutoBuild();
    return cmesh;
}

void ForceFunction(const VecDouble &x, VecDouble &f) {
    f.resize(3);

    double xv = x[0];
    double yv = x[1];
    //    STATE zv = x[2];

    double f_x = 8.0 * Pi * Pi * cos(2.0 * Pi * yv) * sin(2.0 * Pi * xv);
    double f_y = -8.0 * Pi * Pi * cos(2.0 * Pi * xv) * sin(2.0 * Pi * yv);
    double f_z = 0;

    f[0] = f_x; // x direction
    f[1] = f_y; // y direction
    f[2] = f_z; // y direction
}

void Sol_exact(const VecDouble &x, VecDouble &sol, Matrix &dsol) {

    dsol.Resize(3, 3);
    sol.resize(3);

    double xv = x[0];
    double yv = x[1];

    double v_x = cos(2 * Pi * yv) * sin(2 * Pi * xv);
    double v_y = -(cos(2 * Pi * xv) * sin(2 * Pi * yv));
    double v_z = 0;


    sol[0] = v_x;
    sol[1] = v_y;
    sol[2] = v_z;

    // vx direction
    dsol(0, 0) = 2 * Pi * cos(2 * Pi * xv) * cos(2 * Pi * yv);
    dsol(0, 1) = 2 * Pi * sin(2 * Pi * xv) * sin(2 * Pi * yv);
    dsol(0, 2) = 0;

    // vy direction
    dsol(1, 0) = -2 * Pi * sin(2 * Pi * xv) * sin(2 * Pi * yv);
    dsol(1, 1) = -2 * Pi * cos(2 * Pi * xv) * cos(2 * Pi * yv);
    dsol(1, 2) = 0;

    dsol(2, 0) = 0;
    dsol(2, 1) = 0;
    dsol(2, 2) = 0;
}

void QuadrilateralTest(int pOrder) {
    VecDouble error0(3, 0);
    string filename;
    ofstream fout;
    if (pOrder == 1) {
        fout.open("Quadrilateral-Linear.txt");
        filename = "Quadrilateral-Linear.vtk";
    } else if (pOrder == 2) {
        fout.open("Quadrilateral-Quadratic.txt");
        filename = "Quadrilateral-Quadratic.vtk";
    } else {
        cout << "Not implemented for this order!" << endl;
        DebugStop();
    }

    for (int i = 0; i < 5; i++) {
        int ndiv = pow(2, i);
        double l = 1;

        fout << "-------------------------------------------" << std::endl;
        fout << "Number of elements: " << ndiv << "x" << ndiv << std::endl;
        cout << "\nNumber of elements: " << ndiv << "x" << ndiv << std::endl;

        GeoMesh *gmesh = QuadGeoMesh(ndiv + 1, ndiv + 1, l);
        CompMesh *cmesh = CreateCompMesh(gmesh, pOrder, 2);

        Analysis an(cmesh);
        PostProcess *pos = new PostProcessTemplate<Poisson>(&an);
        pos->AppendVariable("Sol");
        pos->AppendVariable("SolExact");
        pos->SetExact(Sol_exact);

        an.RunSimulation();
        an.PostProcessSolution(filename, *pos);

        VecDouble error;
        error = an.PostProcessError(fout, *pos);

        ComputeConvergenceRates(fout, error0, error, ndiv, l);
        fout << "-------------------------------------------" << std::endl;
        delete pos;
        delete gmesh;
        delete cmesh;
        error.clear();
    }
    error0.clear();
}

void TriangleTest(int pOrder) {
    VecDouble error0(3, 0);
    string filename;
    ofstream fout;
    if (pOrder == 1) {
        fout.open("Triangle-Linear.txt");
        filename = "Triangle-Linear.vtk";
    } else if (pOrder == 2) {
        fout.open("Triangle-Quadratic.txt");
        filename = "Triangle-Quadratic.vtk";
    } else {
        cout << "Not implemented for this order!" << endl;
        DebugStop();
    }

    for (int i = 0; i < 5; i++) {
        int ndiv = pow(2, i);
        double l = 1;

        fout << "-------------------------------------------" << std::endl;
        fout << "Number of elements: " << ndiv << "x" << ndiv << std::endl;
        cout << "\nNumber of elements: " << ndiv << "x" << ndiv << std::endl;

        GeoMesh *gmesh = TriangleGeoMesh(ndiv + 1, ndiv + 1, l);
        CompMesh *cmesh = CreateCompMesh(gmesh, pOrder, 2);

        Analysis an(cmesh);
        PostProcess *pos = new PostProcessTemplate<Poisson>(&an);
        pos->AppendVariable("Sol");
        pos->AppendVariable("SolExact");
        pos->SetExact(Sol_exact);

        an.RunSimulation();
        an.PostProcessSolution(filename, *pos);

        VecDouble error;
        error = an.PostProcessError(fout, *pos);

        ComputeConvergenceRates(fout, error0, error, ndiv, l);
        fout << "-------------------------------------------" << std::endl;
        delete pos;
        delete gmesh;
        delete cmesh;
        error.clear();
    }
    error0.clear();
}

void TetrahedronTest(int pOrder) {
    VecDouble error0(3, 0);
    ofstream fout;
    string filename;
    if (pOrder == 1) {
        fout.open("Tetrahedron-Linear.txt");
        filename = "Tetrahedron-Linear.vtk";
    } else if (pOrder == 2) {
        fout.open("Tetrahedron-Quadratic.txt");
        filename = "Tetrahedron-Quadratic.vtk";
    } else {
        cout << "Not implemented for this order!" << endl;
        DebugStop();
    }

    for (int i = 0; i < 1; i++) {
        int ndiv = pow(2, 0);
        double l = 1;

        fout << "-------------------------------------------" << std::endl;
        fout << "Number of elements: " << ndiv << "x" << ndiv << std::endl;
        cout << "\nNumber of elements: " << ndiv << "x" << ndiv << std::endl;

        GeoMesh *gmesh = TetrahedronGeoMesh(ndiv + 1, ndiv + 1, ndiv + 1, l);
        CompMesh *cmesh = CreateCompMesh(gmesh, pOrder, 3);

        Analysis an(cmesh);
        PostProcess *pos = new PostProcessTemplate<Poisson>(&an);
        pos->AppendVariable("Sol");
        pos->AppendVariable("SolExact");
        pos->SetExact(Sol_exact);

        an.RunSimulation();
        an.PostProcessSolution(filename, *pos);

        VecDouble error;
        error = an.PostProcessError(fout, *pos);

        ComputeConvergenceRates(fout, error0, error, ndiv, l);
        fout << "-------------------------------------------" << std::endl;
        delete pos;
        delete gmesh;
        delete cmesh;
        error.clear();
    }
    error0.clear();
}


