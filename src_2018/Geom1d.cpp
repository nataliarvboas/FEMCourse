/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi) {
    double qsi = xi[0];
    
    phi[0] = (1.0 - qsi) / 2.;
    phi[1] = (1.0 + qsi) / 2.;
    
    dphi(0, 0) = -0.5;
    dphi(0, 1) = 0.5;
}

void Geom1d::X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x) {
    int nrow = NodeCo.Rows();
    for (int i = 0; i < nrow; i++) {
        x[i] = NodeCo.GetVal(i, 0)*(1. - xi[0])*0.5 + NodeCo.GetVal(i, 1)*(1. + xi[0])*0.5;
    }
}

void Geom1d::GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx) {
    int nrow = NodeCo.Rows();
    int ncol = NodeCo.Cols();

    gradx.Resize(nrow, 1);
    gradx.Zero();

    VecDouble phi(2);
    Matrix dphi(2, 2);
    Shape(xi, phi, dphi);
    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nrow; j++) {
            gradx(j, 0) += NodeCo.GetVal(j, i) * dphi(0, i);
        }
    }
    phi.clear();
}

void Geom1d::SetNodes(const VecInt &nodes) {
    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt &nodes) {
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) {
    return fNodeIndices[node];
}

int Geom1d::NumNodes() {
    return nCorners;    
}

GeoElementSide Geom1d::Neighbour(int side) {
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side]=neighbour;
}