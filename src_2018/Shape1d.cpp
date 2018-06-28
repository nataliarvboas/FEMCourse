/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Shape1d.h"
//#include "TMatrix.h"

void Shape1d::Shape(const VecDouble& xi, VecInt& orders, VecDouble& phi, Matrix& dphi) {
    int n = NShapeFunctions(orders);
    
    if (orders[0] == 2) {
        phi.resize(n);
        dphi.Resize(2, n);
        n = n - 1;
    }

    for (int i = 0; i < n; i++) {
        phi[i] = 1.;
        dphi(0, i) = 0;
    }

    for (int i = 0; i < n; i++) {
        double epsi = -1. + i * 2. / (n - 1);

        for (int j = 0; j < n; j++) {

            Matrix axdphi(1, n);
            if (i != j) {
                double epsj = -1. + j * 2. / (n - 1);

                phi[i] *= (xi[0] - epsj) / (epsi - epsj);

                axdphi(0, i) = 1 / (epsi - epsj);

                for (int k = 0; k < n; k++) {
                    if (k != i && k != j) {
                        epsj = -1. + k * 2. / (n - 1);
                        axdphi(0, i) *= (xi[0] - epsj) / (epsi - epsj);
                    }
                }
                dphi(0, i) += axdphi(0, i);
            }
        }
    }
    if (orders[0] == 2) {
        phi[2] = phi[0] * phi[1];
        dphi(0, 2) = dphi(0, 0) * phi[1] + phi[0] * dphi(0, 1);
        phi[2] *= 4.;
        dphi(0, 2) *= 4.;
    }

}

int Shape1d::NShapeFunctions(int side, int order) {
    if (side < 2) {
        return 1;
    }
    if (side < 3) {
        return (order - 1);
    }

}

int Shape1d::NShapeFunctions(VecInt &orders) {
    int n = orders.size();
    int nshape = 0;
    int i = 0;
    for (i = 0; i < n; i++) {
        nshape = nshape + NShapeFunctions(i, orders[i]);
    }
    return nshape;
}