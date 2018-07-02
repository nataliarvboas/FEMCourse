/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShapeTetrahedron.h"
#include "tpanic.h"

void ShapeTetrahedron::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {
    int nshape = NShapeFunctions(orders);

    phi[0] = 1 - xi[0] - xi[1] - xi[2];
    phi[1] = xi[0];
    phi[2] = xi[1];
    phi[3] = xi[2];

    dphi(0, 0) = -1.0;
    dphi(0, 1) = 1.0;
    dphi(0, 2) = 0.0;
    dphi(0, 3) = 0.0;

    dphi(1, 0) = -1.0;
    dphi(1, 1) = 0.0;
    dphi(1, 2) = 1.0;
    dphi(1, 3) = 0.0;

    dphi(2, 0) = -1.0;
    dphi(2, 1) = 0.0;
    dphi(2, 2) = 0.0;
    dphi(2, 3) = 1.0;

    if (nshape == 10) {
        int is;
        for (is = 4; is < nshape; is++) {
            int nsnodes = NSideNodes(is);
            switch (nsnodes) {
                case 2:
                {
                    int is1 = SideNodeIndex(is, 0);
                    int is2 = SideNodeIndex(is, 1);
                    phi[is] = phi[is1] * phi[is2];
                    dphi(0, is) = dphi(0, is1) * phi[is2] + phi[is1] * dphi(0, is2);
                    dphi(1, is) = dphi(1, is1) * phi[is2] + phi[is1] * dphi(1, is2);
                    dphi(2, is) = dphi(2, is1) * phi[is2] + phi[is1] * dphi(2, is2);
                }
                    break;
                case 3:
                {
                    //int face = is-10;
                    int is1 = SideNodeIndex(is, 0); //ShapeFaceId[face][0]; 
                    int is2 = SideNodeIndex(is, 1); //ShapeFaceId[face][1]; 
                    int is3 = SideNodeIndex(is, 2); //ShapeFaceId[face][2]; 
                    phi[is] = phi[is1] * phi[is2] * phi[is3];
                    dphi(0, is) = dphi(0, is1) * phi[is2] * phi[is3] + phi[is1] * dphi(0, is2) * phi[is3] + phi[is1] * phi[is2] * dphi(0, is3);
                    dphi(1, is) = dphi(1, is1) * phi[is2] * phi[is3] + phi[is1] * dphi(1, is2) * phi[is3] + phi[is1] * phi[is2] * dphi(1, is3);
                    dphi(2, is) = dphi(2, is1) * phi[is2] * phi[is3] + phi[is1] * dphi(2, is2) * phi[is3] + phi[is1] * phi[is2] * dphi(2, is3);
                }
                    break;
                case 4:
                {
                    phi[is] = phi[0] * phi[1] * phi[2] * phi[3];
                    for (int xj = 0; xj < 3; xj++) {
                        dphi(xj, is) = dphi(xj, 0) * phi[1] * phi[2] * phi[3] +
                                phi[0] * dphi(xj, 1) * phi[2] * phi[3] +
                                phi[0] * phi[1] * dphi(xj, 2) * phi[3] +
                                phi[0] * phi[1] * phi[2] * dphi(xj, 3);
                    }
                }
                    break;

                default:
                    DebugStop();
            }
        }

        double mult[] = {1., 1., 1., 1., 4., 4., 4., 4., 4., 4., 27., 27., 27., 27., 54.};
        for (is = 4; is < nshape; is++) {
            phi[is] *= mult[is];
            dphi(0, is) *= mult[is];
            dphi(1, is) *= mult[is];
            dphi(2, is) *= mult[is];
        }
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