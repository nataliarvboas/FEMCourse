/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Shape1d.h"
#include "ShapeQuad.h"

void ShapeQuad::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {

    int nshape = NShapeFunctions(orders);
    int nsides = orders.size();

    Matrix Indices(2, 2, 0);
    Indices(0, 0) = 0;
    Indices(0, 1) = 3;
    Indices(1, 0) = 1;
    Indices(1, 1) = 2;

    phi.resize(nshape);
    dphi.Resize(2, nshape);

    VecDouble phixi(nshape), phieta(nshape);
    Matrix dphixi(1, nshape), dphieta(1, nshape);

    VecDouble cxi(1);
    cxi[0] = xi[0];
    VecDouble ceta(1);
    ceta[0] = xi[1];

    VecInt order_1d(3, 1);
    for (int i = 0; i < 3; i++) {
        order_1d[i] = 1;
    }

    int nn = Shape1d::NShapeFunctions(order_1d);
    for (int csi = 0; csi < nn; csi++) {

        Shape1d::Shape(cxi, order_1d, phixi, dphixi);

        for (int eta = 0; eta < nn; eta++) {

            Shape1d::Shape(ceta, order_1d, phieta, dphieta);

            phi[Indices(csi, eta)] = phixi[csi] * phieta[eta];

            dphi(0, Indices(csi, eta)) = dphixi(0, csi) * phieta[eta];
            dphi(1, Indices(csi, eta)) = dphieta(0, eta) * phixi[csi];
        }
    }

    if (orders[0] == 2) {
        int is;
        for (is = 4; is < 8; is++) {
            int id = is % 4;
            int id2 = (is + 1) % 4;
            phi[is] = phi[id] * phi[id2];
            dphi(0, is) = dphi(0, id) * phi[id2] + phi[id] * dphi(0, id2);
            dphi(1, is) = dphi(1, id) * phi[id2] + phi[id] * dphi(1, id2);
        }
        phi[8] = phi[0] * phi[2];
        dphi(0, 8) = dphi(0, 0) * phi[2] + phi[0] * dphi(0, 2);
        dphi(1, 8) = dphi(1, 0) * phi[2] + phi[0] * dphi(1, 2);

        for (int is = 4; is < 8; is++) {
            phi[is] += phi[8];
            dphi(0, is) += dphi(0, 8);
            dphi(1, is) += dphi(1, 8);
            phi[is] *= 4.;
            dphi(0, is) *= 4.;
            dphi(1, is) *= 4.;
        }
        phi[8] *= 16.;
        dphi(0, 8) *= 16.;
        dphi(1, 8) *= 16.;
    }
    phixi.clear();
    phieta.clear();
    cxi.clear();
    ceta.clear();
    order_1d.clear();
}

int ShapeQuad::NShapeFunctions(int side, int order) {
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
