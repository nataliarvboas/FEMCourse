/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShapeTetrahedron.h"
#include "TMatrix.h"

void ShapeTetrahedron::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {
    int nshape = NShapeFunctions(orders);

    if (nshape == 4) {
        phi[0] = 1 - xi[0] - xi[1] - xi[2];
        phi[1] = xi[0];
        phi[2] = xi[1];
        phi[3] = xi[2];

        dphi(0, 0) = -1;
        dphi(0, 1) = 1;
        dphi(0, 2) = 0;
        dphi(0, 3) = 0;

        dphi(1, 0) = -1;
        dphi(1, 1) = 0;
        dphi(1, 2) = 1;
        dphi(1, 2) = 0;

        dphi(1, 0) = -1;
        dphi(1, 1) = 0;
        dphi(1, 2) = 0;
        dphi(1, 2) = 1;
    }

    if (nshape == 10) {

        VecDouble eps(3);
        eps[0] = 1 - xi[0] - xi[1] - xi[2];
        eps[1] = xi[0];
        eps[2] = xi[1];
        eps[3] = xi[2];


        for (int i = 0; i <= 3; i++) {
            phi[i] = eps[i]*(2. * eps[i] - 1.);
        }
        phi[4] = 4. * phi[0] * phi[1];
        phi[5] = 4. * phi[1] * phi[2];
        phi[6] = 4. * phi[2] * phi[0];
        phi[7] = 4. * phi[0] * phi[3];
        phi[8] = 4. * phi[1] * phi[3];
        phi[9] = 4. * phi[2] * phi[3];

        dphi(0, 0) = 1. - 4. * (1 - eps[0] - eps[1] - eps[2]);
        dphi(0, 1) = -1. + 4. * eps[0];
        dphi(0, 2) = 0.;
        dphi(0, 3) = 0.;
        dphi(0, 4) = -4. * eps[0]*(-1. + 2. * eps[0])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) - 8. * eps[0]*(-1. + 2. * eps[0])*(1. - eps[0] - eps[1] - eps[2]) + 8. * eps[0]*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2]))*(1. - eps[0] - eps[1] - eps[2]) + 4. * (-1. + 2. * eps[0])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2]))*(1. - eps[0] - eps[1] - eps[2]);
        dphi(0, 5) = 8. * eps[0] * eps[1]*(-1. + 2. * eps[1]) + 4. * (-1. + 2. * eps[0]) * eps[1]*(-1. + 2. * eps[1]);
        dphi(0, 6) = -4 * eps[1]*(-1. + 2. * eps[1])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) - 8. * eps[1]*(-1. + 2. * eps[1])*(1. - eps[0] - eps[1] - eps[2]);
        dphi(0, 7) = -4 * (-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) * eps[2]*(-1. + 2. * eps[2]) - 8. * (1. - eps[0] - eps[1] - eps[2]) * eps[2]*(-1. + 2. * eps[2]);
        dphi(0, 8) = 8. * eps[0] * eps[2]*(-1. + 2. * eps[2]) + 4. * (-1. + 2. * eps[0]) * eps[2]*(-1. + 2. * eps[2]);
        dphi(0, 9) = 0.;

        dphi(1, 0) = 1. - 4. * (1. - eps[0] - eps[1] - eps[2]);
        dphi(1, 1) = 0.;
        dphi(1, 2) = -1. + 4. * eps[1];
        dphi(1, 3) = 0.;
        dphi(1, 4) = -4. * eps[0]*(-1. + 2. * eps[0])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) - 8. * eps[0]*(-1. + 2. * eps[0])*(1. - eps[0] - eps[1] - eps[2]);
        dphi(1, 5) = 8. * eps[0]*(-1. + 2. * eps[0]) * eps[1] + 4 * eps[0]*(-1. + 2. * eps[0])*(-1 + 2. * eps[1]);
        dphi(1, 6) = -4. * eps[1]*(-1. + 2. * eps[1])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) - 8. * eps[1]*(-1. + 2. * eps[1])*(1. - eps[0] - eps[1] - eps[2]) + 8. * eps[1]*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2]))*(1. - eps[0] - eps[1] - eps[2]) + 4. * (-1. + 2. * eps[1])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2]))*(1. - eps[0] - eps[1] - eps[2]);
        dphi(1, 7) = -4. * (-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) * eps[2]*(-1. + 2. * eps[2]) - 8. * (1. - eps[0] - eps[1] - eps[2]) * eps[2]*(-1. + 2. * eps[2]);
        dphi(1, 8) = 0.;
        dphi(1, 9) = 8. * eps[1] * eps[2]*(-1. + 2. * eps[2]) + 4. * (-1. + 2. * eps[1]) * eps[2]*(-1. + 2. * eps[2]);

        dphi(2, 0) = 1. - 4. * (1. - eps[0] - eps[1] - eps[2]);
        dphi(2, 1) = 0.;
        dphi(2, 2) = 0.;
        dphi(2, 3) = -1. + 4. * eps[2];
        dphi(2, 4) = -4. * eps[0]*(-1. + 2. * eps[0])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) - 8. * eps[0]*(-1. + 2. * eps[0])*(1. - eps[0] - eps[1] - eps[2]);
        dphi(2, 5) = 0.;
        dphi(2, 6) = -4 * eps[1]*(-1. + 2. * eps[1])*(-1. + 2. * (1. - eps[0] - eps[1] - eps[2])) - 8. * eps[1]*(-1. + 2. * eps[1])*(1. - eps[0] - eps[1] - eps[2]);
        dphi(2, 7) = 8. * (-1. + 2. * (1. - eps[0] - eps[1] - eps[2]))*(1. - eps[0] - eps[1] - eps[2]) * eps[2] + 4. * (-1. + 2. * (1. - eps[0] - eps[1] - eps[2]))*(1. - eps[0] - eps[1] - eps[2])*(-1. + 2. * eps[2]) - 4. * (-1 + 2. * (1. - eps[0] - eps[1] - eps[2])) * eps[2]*(-1. + 2. * eps[2]) - 8. * (1. - eps[0] - eps[1] - eps[2]) * eps[2]*(-1. + 2. * eps[2]);
        dphi(2, 8) = 8. * eps[0]*(-1. + 2. * eps[0]) * eps[2] + 4. * eps[0]*(-1. + 2. * eps[0])*(-1. + 2. * eps[2]);
        dphi(2, 9) = 8. * eps[1]*(-1. + 2. * eps[1]) * eps[2] + 4. * eps[1]*(-1. + 2. * eps[1])*(-1. + 2. * eps[2]);
    }

}

int ShapeTetrahedron::NShapeFunctions(int side, int order) {
    if (side < 4) return 1;
    if (side < 10) return order - 1;
    if (side < 14) {
        int sum = 0;
        for (int i = 0; i < order - 1; i++) sum += i;
        return sum;
    }
    if (side == 14) {
        int totsum = 0, sum;
        for (int i = 1; i < order - 2; i++) {
            sum = i * (i + 1) / 2;
            totsum += sum;
        }
        return totsum;
    }
}

int ShapeTetrahedron::NShapeFunctions(VecInt &orders) {
    int n = orders.size();
    int nshape = 0;
    for (int i = 0; i < n; i++) {
        nshape += NShapeFunctions(i, orders[i]);
    }
    return nshape;
}