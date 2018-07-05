/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "tpanic.h"

IntRuleQuad::IntRuleQuad(){
}

IntRuleQuad::IntRuleQuad(int order) {
    if (order < 0) {
        DebugStop();
    }

    SetOrder(order);
}

void IntRuleQuad::SetOrder(int order) {
    fOrder = order;
    
        int npoints;
    int resto = fOrder % 2;

    if (resto != 0) {
        npoints = (fOrder + 1) / 2;
    }
    if (resto == 0) {
        npoints = (fOrder + 2) / 2;
    }

    fPoints.Resize(npoints*npoints, 2);
    fWeights.resize(npoints * npoints);
    
    if (fOrder <= 19) {
        IntRule1d x(fOrder);
        IntRule1d y(fOrder);

        double weight;
        VecDouble co(1);

        for (int i = 0; i < npoints; i++) {

            x.Point(i, co, weight);

            VecDouble coX(1);
            double weightX;
            coX[0] = co[0];
            weightX = weight;

            for (int j = 0; j < npoints; j++) {

                y.Point(j, co, weight);

                fPoints(j + i * npoints, 0) = co[0];
                fPoints(j + i * npoints, 1) = coX[0];

                fWeights[j + i * npoints] = weightX*weight;
            }
        }
        co.clear();
    }

    if (order > 19) {
        VecDouble co(npoints);
        VecDouble w(npoints);

        this->gaulegQuad(-1, 1, co, w);

        for (int i = 0; i < npoints; i++) {
            for (int j = 0; j < npoints; j++) {
            fPoints(j + i * npoints, 0) = co[j + i * npoints];
            fPoints(j + i * npoints, 1) = co[j + i * npoints + npoints * npoints];

            fWeights[j + i * npoints] = w[i];
            }
        }
        co.clear();
        w.clear();        
    }   
}

void IntRuleQuad::gaulegQuad(const double x1, const double x2, VecDouble &co, VecDouble &w) {
    IntRule1d x;
    IntRule1d y;
    
    int n = w.size();   

    VecDouble cox(n);
    VecDouble coy(n);
    VecDouble wx(n);
    VecDouble wy(n);


    x.gauleg(x1, x2, cox, wx);
    y.gauleg(x1, x2, coy, wy);
    
    co.resize(2*n*n);
    w.resize(n * n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            co[j + i * n] = cox[j];
            co[j + i * n + n * n] = coy[i];
            w[n * i + j] = wx[i] * wy[j];
        }
    }
    cox.clear();
    coy.clear();
    wx.clear();
    wy.clear();    
}