/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "TMatrix.h"

void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {

    int nshape = NShapeFunctions(orders);
    int nsides = orders.size();

    Matrix Indices(1, 1, 0);

    if (nsides == 4) {

        Indices.Resize(2, 2);

        Indices(0, 0) = 0;
        Indices(0, 1) = 3;
        Indices(1, 0) = 1;
        Indices(1, 1) = 2;

        phi.resize(4);
        dphi.Resize(2, 4);
    }

    if (nsides == 9) {
        Indices.Resize(3, 3);

        Indices(0, 0) = 0;
        Indices(0, 1) = 7;
        Indices(0, 2) = 3;
        Indices(1, 0) = 4;
        Indices(1, 1) = 8;
        Indices(1, 2) = 6;
        Indices(2, 0) = 1;
        Indices(2, 1) = 5;
        Indices(2, 2) = 2;

        phi.resize(9);
        dphi.Resize(2, 9);
    }

    VecDouble coxi(1);
    coxi[0] = xi[0];

    VecDouble coeta(1);
    coeta[0] = xi[1];

    VecDouble phixi(nshape), phieta(nshape);
    TMatrix dphixi(1, nshape), dphieta(1, nshape);

    for (int csi = 0; csi < nshape; csi++) {

        Shape1d::Shape(xi, orders, phixi, dphixi);

        for (int eta = 0; eta < nshape; eta++) {

            Shape1d::Shape(xi, orders, phixi, dphixi);

            phi[Indices(csi, eta)] = phixi[csi] * phieta[eta];

            dphi(0, Indices(csi, eta)) = dphixi(0, csi) * phieta[eta];
            dphi(1, Indices(csi, eta)) = dphieta(0, eta) * phixi[csi];
        }
    }
}

int ShapeQuad::NShapeFunctions(int side, int order) {
    if (side < 4) return 1;
    return order - 1;

    if (side < 4) return 1;

    if (side < 8) return (order - 1); 
    if (side == 8) {
        return ((order - 1)*(order - 1));
    }

}

int ShapeQuad::NShapeFunctions(VecInt &orders) {
    int n = orders.size();
    int nshape = 0;
    for (int i = 0; i < n; i++) {
        nshape += NShapeFunctions(i, orders[i]);
    }
    return nshape;
}
