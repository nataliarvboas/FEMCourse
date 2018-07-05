/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeomTetrahedron.h"

GeomTetrahedron::GeomTetrahedron() {
}

GeomTetrahedron::~GeomTetrahedron() {
}

GeomTetrahedron::GeomTetrahedron(const GeomTetrahedron &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomTetrahedron& GeomTetrahedron::operator=(const GeomTetrahedron& copy){
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomTetrahedron::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    double qsi = xi[0];
    double eta = xi[1];
    double zeta = xi[2];

    phi[0] = 1.0 - qsi - eta - zeta;
    phi[1] = qsi;
    phi[2] = eta;
    phi[3] = zeta;

    dphi(0, 0) = -1.0;
    dphi(1, 0) = -1.0;
    dphi(2, 0) = -1.0;
    
    dphi(0, 1) = 1.0;
    dphi(1, 1) = 0.0;
    dphi(2, 1) = 0.0;
   
    dphi(0, 2) = 0.0;
    dphi(1, 2) = 1.0;
    dphi(2, 2) = 0.0;
    
    dphi(0, 3) = 0.0;
    dphi(1, 3) = 0.0;
    dphi(2, 3) = 1.0;
}

void GeomTetrahedron::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    VecDouble phi(4);
    Matrix dphi(3, 4);
    
    
    Shape(xi, phi, dphi);
    int space = NodeCo.Rows();

    for (int i = 0; i < space; i++) {
        x[i] = 0.0;
        for (int j = 0; j < 4; j++) {
            x[i] += phi[j] * NodeCo.GetVal(i, j);
        }
    }
    phi.clear();
}


void GeomTetrahedron::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    gradx.Resize(3, 3);
    gradx.Zero();
    int nrow = NodeCo.Rows();
    int ncol = NodeCo.Cols();

    VecDouble phi(4);
    Matrix dphi(3, 4);
    Shape(xi, phi, dphi);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            gradx(j, 0) += NodeCo.GetVal(j, i) * dphi(0, i);
            gradx(j, 1) += NodeCo.GetVal(j, i) * dphi(1, i);
            gradx(j, 2) += NodeCo.GetVal(j, i) * dphi(2, i);
        }
    }
    phi.clear();
}


void GeomTetrahedron::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void GeomTetrahedron::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int GeomTetrahedron::NodeIndex(int node) {
    return fNodeIndices[node];
}

int GeomTetrahedron::NumNodes() {
    return nCorners;

}

GeoElementSide GeomTetrahedron::Neighbour(int side) {
    return fNeighbours[side];
}

void GeomTetrahedron::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side]=neighbour;
}