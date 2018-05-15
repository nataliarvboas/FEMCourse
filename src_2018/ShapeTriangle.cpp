/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShapeTriangle.h"
#include "TMatrix.h"

void ShapeTriangle::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {

    int nshape = NShapeFunctions(orders);
    int nsides = orders.size();

    phi.resize(nshape);
    dphi.Resize(2, nshape);

    if (nsides == 3) {
        phi[0] = 1 - xi[0] - xi[1];
        phi[1] = xi[0];
        phi[2] = xi[1];


        dphi(0, 0) = -1;
        dphi(0, 1) = 1;
        dphi(0, 2) = 0;
        dphi(1, 0) = -1;
        dphi(1, 1) = 0;
        dphi(1, 2) = 1;
    }

    if (nsides == 7) {
        VecDouble eps(nshape);

        eps[0] = 1 - xi[0] - xi[1];
        eps[1] = xi[0];
        eps[2] = xi[1];

        phi[0] = 2 * eps[0]*(eps[0] - 0.5);
        phi[1] = 2 * eps[1]*(eps[1] - 0.5);
        phi[2] = 2 * eps[2]*(eps[2] - 0.5);
        phi[3] = 4 * eps[0] * eps[1];
        phi[4] = 4 * eps[1] * eps[2];
        phi[5] = 4 * eps[2] * eps[0];

        dphi(0, 0) = -2 * (0.5 - xi[1] - xi[0]) - 2 * (1 - xi[1] - xi[0]);
        dphi(0, 1) = 2 * (-0.5 + xi[0]) + 2 * xi[0];
        dphi(0, 2) = 0;
        dphi(0, 3) = 4 * (1 - xi[1] - xi[0]) - 4 * xi[0];
        dphi(0, 4) = 4 * xi[1];
        dphi(0, 5) = -4 * xi[1];

        dphi(1, 0) = -2 * (0.5 - xi[1] - xi[0]) - 2 * (1 - xi[1] - xi[0]);
        dphi(1, 1) = 0;
        dphi(1, 2) = 2 * (-0.5 + xi[1]) + 2 * xi[1];
        dphi(1, 3) = -4 * xi[0];
        dphi(1, 4) = 4 * xi[0];
        dphi(1, 5) = -4 * xi[1] + 4 * (1 - xi[1] - xi[0]);
    }
}

int ShapeTriangle::NShapeFunctions(int side, int order) {
    if (side < 3) return 1;
    return order - 1;
    switch (side) {
        case 0:
        case 1:
        case 2:
            return 1;
        case 3:
        case 4:
        case 5:
            return order - 1;
        case 6:
            return (order - 2) < 0 ? 0 : ((order - 2)*(order - 1)) / 2;
    }
}

int ShapeTriangle::NShapeFunctions(VecInt &orders) {
    int n = orders.size();
    int nshape = 0;
    for (int i = 0; i < n; i++) {
        nshape += NShapeFunctions(i, orders[i]);
    }
    return nshape;
}
