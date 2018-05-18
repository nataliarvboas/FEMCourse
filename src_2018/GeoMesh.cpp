/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include "GeoMesh.h"
#include "GeoNode.h"
#include "vector"
#include "GeoElementSide.h"
#include <vector>
#include "tpanic.h"

GeoMesh::GeoMesh() {
}

GeoMesh::GeoMesh(const GeoMesh &copy) {
    Elements = copy.Elements;
    Nodes = copy.Nodes;
}

GeoMesh &GeoMesh::operator=(const GeoMesh &copy) {

    int nnodes = copy.Nodes.size();
    this->Nodes.resize(nnodes);
    for (int64_t inodes = 0; inodes < nnodes; inodes++)
        this->Nodes[inodes] = copy.Nodes[inodes];

    this->Elements.resize(copy.Elements.size());
    for (int64_t iel = 0; iel < copy.Elements.size(); iel++)
        this->Elements[iel] = copy.Elements[iel];

    return *this;
}

void GeoMesh::SetNumNodes(int nnodes) {
    Nodes.resize(nnodes);
}

void GeoMesh::SetNumElements(int numelements) {
    Elements.resize(numelements);
}

int GeoMesh::NumNodes() {
    return Nodes.size();
}

int GeoMesh::NumElements() {
    return Elements.size();
}

GeoNode &GeoMesh::Node(int node) {
    return Nodes[node];
}

void GeoMesh::SetElement(int elindex, GeoElement *gel) {
    Elements[elindex] = gel;
}

GeoElement *GeoMesh::Element(int elindex) {
    return Elements[elindex];
}

void GeoMesh::BuildConnectivity() {
    VecInt sides(this->NumNodes(), -1);
    VecInt vetor(this->NumNodes(), -1);


    int64_t nelem = this->NumElements();
    int64_t iel = 0;

    for (iel = 0; iel < nelem; iel++) {
        GeoElement *gel = Elements[iel];
        if (!gel) continue;
        int ncor = gel->NCornerNodes();
        int in = 0;
        for (in = 0; in < ncor; in++) {
            int64_t nodeindex = gel->NodeIndex(in);
            if (vetor[nodeindex] == -1) {
                vetor[nodeindex] = iel;
                sides[nodeindex] = in;

            } else {
                GeoElementSide one(gel, in);
                GeoElementSide two(Element(vetor[nodeindex]), sides[nodeindex]);

                GeoElementSide neighbour = one.Neighbour();
                if (neighbour.Element() == 0) DebugStop();

                if (!two.IsNeighbour(one)) {
                    one.InsertConnectivity(two);
                }
            }
        }
    }

    for (iel = 0; iel < nelem; iel++) {
        GeoElement *gel = Elements[iel];
        if (!gel) continue;
        int ncor = gel->NCornerNodes();
        int nsides = gel->NSides();
        int is;
        for (is = ncor; is < nsides; is++) {
                GeoElementSide gelside(gel, is);
                std::vector<GeoElementSide> neighbours;
                gelside.ComputeNeighbours(neighbours);
                int64_t nneigh = neighbours.size();
                int64_t in;
                for (in = 0; in < nneigh; in++) {
                    gelside.InsertConnectivity(neighbours[in]);
                }
        }
    }
}

void GeoMesh::Print(std::ostream &out) {
    out << "\n------------Mesh Information------------" << std::endl;
    out << "Number of nodes: " << this->NumNodes() << std::endl;
    out << "Number of elements: " << this->NumElements() << std::endl;

    out << "\n-------Geometric Node Information-------" << std::endl;
    for (int i = 0; i < this->NumNodes(); i++) {
        out << "Node index: " << i << "\t\t";
        this->Node(i).Print(out);
    }

    out << "\n------Geometric Element Information------" << std::endl;
    for (int i = 0; i < this->NumElements(); i++) {
        out << "Element index:\t\t" << i << std::endl;
        this->Element(i)->Print(out);
        out << std::endl;
    }
}