
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

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

void ComputeConvergenceRates(VecDouble &error0, VecDouble &error, int ndiv, double l);
GeoMesh *CreateGeoMesh(int nel_x, int nel_y, int dim, double l_x, double l_y);
CompMesh *CreateCompMesh(GeoMesh *gmesh, int pOrder);
void ForceFunction(const VecDouble &co, VecDouble &result);
void Sol_exact(const VecDouble &x, VecDouble &sol, Matrix &dsol);

int main() {
    VecDouble error0(3,0);
    for (int i = 1; i < 5; i++) {

        int ndiv = pow(2, i);
        int dim = 2;
        int nel_x = ndiv;
        int nel_y = ndiv;
        double l = 1;
        double l_x = l;
        double l_y = l;
        int pOrder = 1;

        GeoMesh *gmesh = CreateGeoMesh(nel_x, nel_y, dim, l_x, l_y);
        CompMesh *cmesh = CreateCompMesh(gmesh, pOrder);

        Analysis an(cmesh);
        PostProcess *pos = new PostProcessTemplate<Poisson>(&an);
        pos->AppendVariable("Sol");
        pos->AppendVariable("SolExact");
        pos->SetExact(Sol_exact);

        an.RunSimulation();
        an.PostProcessSolution("sol.vtk", *pos);
        
        VecDouble error;
        error = an.PostProcessError(std::cout, *pos);
   
        ComputeConvergenceRates(error0, error, ndiv, l);
    }
    return 0;
}

void ComputeConvergenceRates(VecDouble &error0, VecDouble &error, int ndiv, double l) {

    if (ndiv == 2) {
        error0 = error;
    }

    double h = l / ndiv;
    double h0 = 2 * h;

    std::cout << "-----------------------" << std::endl;
    std::cout << "Taxa de conv. u: " << (log(error[0]) - log(error0[0])) / (log(h) - log(h0)) << endl;
    std::cout << "Taxa de conv. gradu: " << (log(error[1]) - log(error0[1])) / (log(h) - log(h0)) << endl;
    std::cout << "Taxa de conv. energia: " << (log(error[2]) - log(error0[2])) / (log(h) - log(h0)) << endl;
    std::cout << "-----------------------" << std::endl;

    for (int i = 0; i < error.size(); i++) {
        error0[i] = error[i];
    }
}

GeoMesh *CreateGeoMesh(int nel_x, int nel_y, int dim, double l_x, double l_y) {
    int nnodes_x = nel_x + 1;
    int nnodes_y = nel_y + 1;
    int64_t nelem = nel_x*nel_y;
    int64_t nnodes = nnodes_x*nnodes_y;

    GeoMesh *gmesh = new GeoMesh;

    gmesh->SetDimension(dim);
    gmesh->SetNumNodes(nnodes);

    //Cria os nos da malha
    VecDouble coord(3, 0.);
    int64_t nodeid;
    for (int i = 0; i < nnodes_y; i++) {
        for (int j = 0; j < nnodes_x; j++) {
            nodeid = i * nnodes_x + j;
            coord[0] = (j) * l_x / (nnodes_x - 1);
            coord[1] = (i) * l_y / (nnodes_y - 1);
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

            int matid = 0;
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
            for (int j = 0; j < 4; j++) {
                co(i, j) = gmesh->Node(TopolQuad[j]).Co()[i];
            }
        }

        // Condicao de contorno inferior -> sides: 0 e 1
        if (co(1, 0) == 0 && co(1, 1) == 0) {
            topolLine[0] = TopolQuad[0];
            topolLine[1] = TopolQuad[1];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc2, gmesh, id);
            id++;
        }
        // Condicao de contorno da direita -> sides: 1 e 2
        if (co(0, 1) == l_x && co(0, 2) == l_x) {
            topolLine[0] = TopolQuad[1];
            topolLine[1] = TopolQuad[2];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc1, gmesh, id);
            id++;
        }
        // Condicao de contorno superior -> sides: 2 e 3
        if (co(1, 2) == l_y && co(1, 3) == l_y) {
            topolLine[0] = TopolQuad[2];
            topolLine[1] = TopolQuad[3];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc3, gmesh, id);
            id++;
        }

        // Condicao de contorno da esquerda -> sides: 0 e 3
        if (co(0, 0) == 0 && co(0, 3) == 0) {
            topolLine[0] = TopolQuad[0];
            topolLine[1] = TopolQuad[3];
            GeoElement *gel = new GeoElementTemplate<Geom1d> (topolLine, bc0, gmesh, id);
            id++;
        }
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

CompMesh *CreateCompMesh(GeoMesh *gmesh, int pOrder) {
    CompMesh *cmesh = new CompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    int nelem = cmesh->GetElementVec().size();

    Matrix perm(2, 2);
    perm(0, 0) = 1;
    perm(0, 1) = 0;
    perm(1, 0) = 0;
    perm(1, 1) = 1;

    for (int i = 0; i < nelem; i++) {
        int matid = gmesh->Element(i)->Material();
        if (matid == 0) {
            cmesh->SetNumberMath(i + 1);
            Poisson *p = new Poisson(matid, perm);
            p->SetForceFunction(ForceFunction);
            p->SetExactSolution(Sol_exact);
            cmesh->SetMathStatement(i, p);
        } else {
            cmesh->SetNumberMath(i + 1);
            Matrix Val1(1, 1, 0.), Val2(1, 1, 0.);
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
    f.resize(2);

    double xv = x[0];
    double yv = x[1];
    //    STATE zv = x[2];

    double f_x = 8.0 * Pi * Pi * cos(2.0 * Pi * yv) * sin(2.0 * Pi * xv);
    double f_y = -8.0 * Pi * Pi * cos(2.0 * Pi * xv) * sin(2.0 * Pi * yv);

    f[0] = f_x; // x direction
    f[1] = f_y; // y direction
}

void Sol_exact(const VecDouble &x, VecDouble &sol, Matrix &dsol) {

    dsol.Resize(2, 2);
    sol.resize(2);

    double xv = x[0];
    double yv = x[1];

    double v_x = cos(2 * Pi * yv) * sin(2 * Pi * xv);
    double v_y = -(cos(2 * Pi * xv) * sin(2 * Pi * yv));


    sol[0] = v_x;
    sol[1] = v_y;

    // vx direction
    dsol(0, 0) = 2 * Pi * cos(2 * Pi * xv) * cos(2 * Pi * yv);
    dsol(0, 1) = 2 * Pi * sin(2 * Pi * xv) * sin(2 * Pi * yv);

    // vy direction
    dsol(1, 0) = -2 * Pi * sin(2 * Pi * xv) * sin(2 * Pi * yv);
    dsol(1, 1) = -2 * Pi * cos(2 * Pi * xv) * cos(2 * Pi * yv);

}
