/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include "GeoElementTemplate.h"
#include "GeoNode.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTriangle.h"
#include "GeomTetrahedron.h"
#include "GeoElement.h"
#include "GeoElementSide.h"
#include "GeoMesh.h"
#include "tpanic.h"

template<class TGeom>
GeoElementTemplate<TGeom>::GeoElementTemplate(const VecInt &nodeindices, int materialid, GeoMesh *gmesh, int index): GeoElement(materialid, gmesh, index) {
    Geom.SetNodes(nodeindices);
    for(int side = 0; side < TGeom::nSides; side++){
        Geom.SetNeighbour(side,GeoElementSide(this,side));
    }
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
    GeoMesh *gmesh = new GeoMesh;

    int NNodes = gmesh->NumNodes();

    Matrix coord(3, NNodes);

    GeoNode np;

    int i, j;
    for (i = 0; i < NNodes; i++) {
        np = gmesh->Node(i);
        for (j = 0; j < 3; j++) {
            coord(j, i) = np.Coord(j);
        }
    }
    Geom.X(xi, coord, x);
}

template<class TGeom>
void GeoElementTemplate<TGeom>::GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx) {
    GeoMesh *gmesh = new GeoMesh;

    int NNodes = gmesh->NumNodes();

    Matrix coord(3, NNodes);

    GeoNode np;

    int i, j;
    for (i = 0; i < NNodes; i++) {
        np = gmesh->Node(i);
        for (j = 0; j < 3; j++) {
            coord(j, i) = np.Coord(j);
        }
    }
    Geom.GradX(xi, coord, x, gradx);
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

template class GeoElementTemplate<Geom1d>;
template class GeoElementTemplate<GeomQuad>;
template class GeoElementTemplate<GeomTriangle>;
template class GeoElementTemplate<GeomTetrahedron>;