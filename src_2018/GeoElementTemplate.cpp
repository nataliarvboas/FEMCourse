/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include "GeoElementTemplate.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeoMesh.h"
#include "tpanic.h"
#include <math.h> 
using namespace std;

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index) : GeoElement(materialid, gmesh, index) {
    gmesh->SetNumElements(index + 1);
    Geom.SetNodes(nodeindices);
    for (int side = 0; side < TGeom::nSides; side++) {
        Geom.SetNeighbour(side, GeoElementSide(this, side));
    }
    gmesh->SetElement(index, this);
}

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const GeoElementTemplate &copy) {
    Geom = copy.Geom;
}

template<class TGeom>
GeoElementTemplate<TGeom> &GeoElementTemplate<TGeom>::operator=(const GeoElementTemplate &copy) {
    Geom = copy.Geom;
    return *this;
}

template<class TGeom>
ElementType GeoElementTemplate<TGeom>::Type() {
    return TGeom::Type();
}

template<class TGeom>
void GeoElementTemplate<TGeom>::X(const VecDouble &xi, VecDouble &x) {
    int NNodes = this->NNodes();
    Matrix coord(3, NNodes);

    int i, j;
    for (i = 0; i < NNodes; i++) {
        int index = this->NodeIndex(i);
        GeoNode node = GMesh->Node(index);
        for (j = 0; j < 3; j++) {
            coord(j, i) = node.Coord(j);
        }
    }
    Geom.X(xi, coord, x);
}

template<class TGeom>
void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx) {
    int NNodes = this->NNodes();
    Matrix coord(3, NNodes);

    int i, j;
    for (i = 0; i < NNodes; i++) {
        int index = this->NodeIndex(i);
        GeoNode node = GMesh->Node(index);
        for (j = 0; j < 3; j++) {
            coord(j, i) = node.Coord(j);
        }
    }
    Geom.GradX(xi, coord, x, gradx);
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Jacobian(const Matrix &gradx, Matrix &jac, Matrix &axes, double &detjac, Matrix &jacinv) {
    int nrows = gradx.Rows();
    int ncols = gradx.Cols();
    int dim = jac.Rows();
    jac.Resize(dim, dim);
    axes.Resize(dim, 3);
    jacinv.Resize(dim, dim);
    jac.Zero();

    switch (dim) {
        case 1:
        {
            axes.Resize(dim, 3);
            VecDouble v_1(3, 0.);

            for (int i = 0; i < nrows; i++) {
                v_1[i] = gradx.GetVal(i, 0);
            }

            double norm_v_1 = 0.;
            for (int i = 0; i < nrows; i++) {
                norm_v_1 += v_1[i] * v_1[i];
            }

            norm_v_1 = sqrt(norm_v_1);
            jac(0, 0) = norm_v_1;
            detjac = norm_v_1;
            jacinv(0, 0) = 1.0 / detjac;

            detjac = fabs(detjac);

            for (int i = 0; i < 3; i++) {
                axes(0, i) = v_1[i] / norm_v_1;
            }
            v_1.clear();
        }
            break;
        case 2:
        {
            VecDouble v_1(3, 0.), v_2(3, 0.);
            VecDouble v_1_til(3, 0.), v_2_til(3, 0.);

            for (int i = 0; i < nrows; i++) {
                v_1[i] = gradx.GetVal(i, 0);
                v_2[i] = gradx.GetVal(i, 1);
            }

            double norm_v_1_til = 0.0;
            double norm_v_2_til = 0.0;
            double v_1_dot_v_2 = 0.0;

            for (int i = 0; i < 3; i++) {
                norm_v_1_til += v_1[i] * v_1[i];
                v_1_dot_v_2 += v_1[i] * v_2[i];
            }
            norm_v_1_til = sqrt(norm_v_1_til);

            for (int i = 0; i < 3; i++) {
                v_1_til[i] = v_1[i] / norm_v_1_til;
                v_2_til[i] = v_2[i] - v_1_dot_v_2 * v_1_til[i] / norm_v_1_til;
                norm_v_2_til += v_2_til[i] * v_2_til[i];
            }
            norm_v_2_til = sqrt(norm_v_2_til);


            jac(0, 0) = norm_v_1_til;
            jac(0, 1) = v_1_dot_v_2 / norm_v_1_til;
            jac(1, 1) = norm_v_2_til;

            detjac = jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1);

            jacinv(0, 0) = +jac(1, 1) / detjac;
            jacinv(1, 1) = +jac(0, 0) / detjac;
            jacinv(0, 1) = -jac(0, 1) / detjac;
            jacinv(1, 0) = -jac(1, 0) / detjac;

            detjac = fabs(detjac);

            for (int i = 0; i < 3; i++) {
                v_2_til[i] /= norm_v_2_til;
                axes(0, i) = v_1_til[i];
                axes(1, i) = v_2_til[i];
            }
            v_1.clear();
            v_2.clear();
            v_1_til.clear();
            v_2_til.clear();
        }
            break;
        case 3:
        {
            axes.Resize(dim, 3);

            for (int i = 0; i < nrows; i++) {
                jac(i, 0) = gradx.GetVal(i, 0);
                jac(i, 1) = gradx.GetVal(i, 1);
                jac(i, 2) = gradx.GetVal(i, 2);
            }

            detjac -= jac(0, 2) * jac(1, 1) * jac(2, 0); //- a02 a11 a20
            detjac += jac(0, 1) * jac(1, 2) * jac(2, 0); //+ a01 a12 a20
            detjac += jac(0, 2) * jac(1, 0) * jac(2, 1); //+ a02 a10 a21
            detjac -= jac(0, 0) * jac(1, 2) * jac(2, 1); //- a00 a12 a21
            detjac -= jac(0, 1) * jac(1, 0) * jac(2, 2); //- a01 a10 a22
            detjac += jac(0, 0) * jac(1, 1) * jac(2, 2); //+ a00 a11 a22

            jacinv(0, 0) = (-jac(1, 2) * jac(2, 1) + jac(1, 1) * jac(2, 2)) / detjac; //-a12 a21 + a11 a22
            jacinv(0, 1) = (jac(0, 2) * jac(2, 1) - jac(0, 1) * jac(2, 2)) / detjac; //a02 a21 - a01 a22
            jacinv(0, 2) = (-jac(0, 2) * jac(1, 1) + jac(0, 1) * jac(1, 2)) / detjac; //-a02 a11 + a01 a12
            jacinv(1, 0) = (jac(1, 2) * jac(2, 0) - jac(1, 0) * jac(2, 2)) / detjac; //a12 a20 - a10 a22
            jacinv(1, 1) = (-jac(0, 2) * jac(2, 0) + jac(0, 0) * jac(2, 2)) / detjac; //-a02 a20 + a00 a22
            jacinv(1, 2) = (jac(0, 2) * jac(1, 0) - jac(0, 0) * jac(1, 2)) / detjac; //a02 a10 - a00 a12
            jacinv(2, 0) = (-jac(1, 1) * jac(2, 0) + jac(1, 0) * jac(2, 1)) / detjac; //-a11 a20 + a10 a21
            jacinv(2, 1) = (jac(0, 1) * jac(2, 0) - jac(0, 0) * jac(2, 1)) / detjac; //a01 a20 - a00 a21
            jacinv(2, 2) = (-jac(0, 1) * jac(1, 0) + jac(0, 0) * jac(1, 1)) / detjac; //-a01 a10 + a00 a11

            detjac = fabs(detjac);

            axes.Zero();
            axes(0, 0) = 1.0;
            axes(1, 1) = 1.0;
            axes(2, 2) = 1.0;
        }
            break;
    }

}

template<class TGeom>
int GeoElementTemplate<TGeom>::WhichSide(VecInt &SideNodeIds) {
    int64_t cap = SideNodeIds.size();
    int nums = NSides();
    for (int side = 0; side < nums; side++) {
        if (NSideNodes(side) == 2 && cap == 2) {
            int64_t isn1 = NodeIndex(SideNodeIndex(side, 0));
            int64_t isn2 = NodeIndex(SideNodeIndex(side, 1)); //sao = para side<3
            if ((isn1 == SideNodeIds[0] && isn2 == SideNodeIds[1]) ||
                    (isn2 == SideNodeIds[0] && isn1 == SideNodeIds[1])) return side;
        } else if (NSideNodes(side) == 1 && cap == 1) {
            if (NodeIndex(SideNodeIndex(side, 0)) == SideNodeIds[0]) return side;
            //completar
        } else if (NSideNodes(side) == 3 && cap == 3) {
            int64_t sni[3], snx[3], k;
            for (k = 0; k < 3; k++) snx[k] = NodeIndex(SideNodeIndex(side, k)); //el atual
            for (k = 0; k < 3; k++) sni[k] = SideNodeIds[k]; //el viz
            for (k = 0; k < 3; k++) {
                if (snx[0] == sni[k] && snx[1] == sni[(k + 1) % 3] && snx[2] == sni[(k + 2) % 3]) return side;
                if (snx[0] == sni[k] && snx[1] == sni[(k + 2) % 3] && snx[2] == sni[(k + 1) % 3]) return side;
            }//012 120 201 , 021 102 210
        } else if (NSideNodes(side) == 4 && cap == 4) {//face quadrilateral
            int64_t sni[4], snx[4], k;
            for (k = 0; k < 4; k++) snx[k] = NodeIndex(SideNodeIndex(side, k)); //el atual
            for (k = 0; k < 4; k++) sni[k] = SideNodeIds[k]; //vizinho
            if (snx[0] == sni[0]) {
                for (k = 1; k < 4; k++) {
                    if (snx[1] == sni[k] && snx[2] == sni[k % 3 + 1] && snx[3] == sni[(k + 1) % 3 + 1]) return side;
                    if (snx[1] == sni[k] && snx[2] == sni[(k + 1) % 3 + 1] && snx[3] == sni[k % 3 + 1]) return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */
            } else if (snx[1] == sni[0]) {
                for (k = 1; k < 4; k++) {
                    if (snx[0] == sni[k] && snx[2] == sni[k % 3 + 1] && snx[3] == sni[(k + 1) % 3 + 1]) return side;
                    if (snx[0] == sni[k] && snx[2] == sni[(k + 1) % 3 + 1] && snx[3] == sni[k % 3 + 1]) return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 1023 1230 1302 , 1032 1203 1320 */
            } else if (snx[2] == sni[0]) {
                for (k = 0; k < 4; k++) {
                    if (snx[0] == sni[k] && snx[1] == sni[k % 3 + 1] && snx[3] == sni[(k + 1) % 3 + 1]) return side;
                    if (snx[0] == sni[k] && snx[1] == sni[(k + 1) % 3 + 1] && snx[3] == sni[k % 3 + 1]) return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               /* 2013 2130 2301 , 2031 2103 2310 */
            } else if (snx[3] == sni[0]) {
                for (k = 0; k < 4; k++) {
                    if (snx[0] == sni[k] && snx[1] == sni[k % 3 + 1] && snx[2] == sni[(k + 1) % 3 + 1]) return side;
                    if (snx[0] == sni[k] && snx[1] == sni[(k + 1) % 3 + 1] && snx[2] == sni[k % 3 + 1]) return side;
                }// /* 0123 0231 0312 , 0132 0213 0321 */                               / * 3012 3120 3201 , 3021 3102 3210 * /
            }
        } else if (cap < 1 || cap > 4) {
            int is;
            for (is = 0; is < nums; is++) {
                if (NSideNodes(is) == cap) {
                    break;
                }
            }
            if (is != nums) {
                std::cout << "WhichSide must be extended\n";
                DebugStop();
            }
        }
    }
    return -1;
}

template<class TGeom>
void GeoElementTemplate<TGeom>::Print(std::ostream &out) {
    out << "Number of nodes:\t" << this->NNodes() << std::endl;
    out << "Corner nodes:\t\t" << this->NCornerNodes() << std::endl;
    out << "Nodes indexes\t\t";

    for (int i = 0; i < this->NNodes(); i++) out << this->NodeIndex(i) << "   ";
    out << "\nNumber of sides:\t" << this->NSides() << std::endl;
    out << "Material Id:\t\t" << this->Material() << std::endl;

    int nsides = this->NSides();

    for (int i = 0; i < nsides; i++) {
        out << "Neighbours for side " << i << ":\t";
        GeoElementSide neighbour = Neighbour(i);
        GeoElementSide thisside(this, i);
        if (!(neighbour.Element() != 0 && neighbour.Side() > -1)) {
            out << "No neighbour\n";
        } else {
            while (neighbour.Element() != thisside.Element() || neighbour.Side() != thisside.Side()) {
                out << neighbour.Element()->GetIndex() << "/" << neighbour.Side() << ' ';

                neighbour = neighbour.Neighbour();
            }
            out << std::endl;
        }
    }
}

template<class TGeom>
int GeoElementTemplate<TGeom>::SideIsUndefined(int side) {
    if (side < 0 || side > NSides()) {
        std::cout << "GeoElementTemplate<TGeom>::SideIsUndefined - bad side: " << side << std::endl;
        DebugStop();
    }
    return (Geom.Neighbour(side).Side() == -1);
}


template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<GeomTetrahedron>;