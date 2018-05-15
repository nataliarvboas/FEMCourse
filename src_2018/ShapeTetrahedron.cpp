/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "ShapeTetrahedron.h"
#include "TMatrix.h"

void ShapeTetrahedron::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, Matrix &dphi) {

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