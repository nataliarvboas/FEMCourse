/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRule1d.h"
#include <vector>
#include <math.h>
#include <cmath>
#include "tpanic.h"
using namespace std;

#define PI 3.141592654

IntRule1d::IntRule1d(){

}

IntRule1d::IntRule1d(int order) {
    if (order < 0) {
        DebugStop();
    }

    SetOrder(order);

}

void IntRule1d::gauleg(const double x1, const double x2, VecDouble &co, VecDouble &w){
    int n = w.size();

    double EPS = 1.0e-14;
    int m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;    
    
    m = (n + 1) / 2;
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);
    for (i = 0; i < m; i++) {
        z = cos(PI * (i + 0.75) / (n + 0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (j = 0; j < n; j++) {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > EPS);
        co[i] = xm - xl*z;
        co[n - 1 - i] = xm + xl*z;
        w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w[n - 1 - i] = w[i];
    }
    
}

void IntRule1d::SetOrder(int order) {
    fOrder = order;  
    
        if (order < 0) {
        DebugStop();
    }

    int npoints;
    int resto = fOrder % 2;

    if (resto == 0) {
        npoints = (fOrder + 2) / 2;
    } else {
        npoints = (fOrder + 1) / 2;
    }
    
    fPoints.Resize(npoints, 1);
    fWeights.resize(npoints);
    
    if (fOrder == 0 || fOrder == 1) {
        fPoints(0,0) = 0.0000000000;        fWeights[0] = 2.0000000000;
    }
    if (fOrder == 2 || fOrder == 3) {            
        fPoints(0,0) = -0.5773502692;        fWeights[0] = 1.0000000000;
        fPoints(1,0) = 0.5773502692;        fWeights[1] = 1.0000000000;
    }
    if (fOrder == 4 || fOrder == 5) {
        fPoints(0,0) = -0.7745966692;        fWeights[0] = 0.5555555556;
        fPoints(1,0) = 0.0000000000;        fWeights[1] = 0.8888888889;
        fPoints(2,0) = 0.7745966692;        fWeights[2] = 0.5555555556;
    }
    if (fOrder == 6 || fOrder == 7) {        
        fPoints(0,0) = -0.8611363116;        fWeights[0] = 0.3478548451;
        fPoints(1,0) = -0.3399810436;        fWeights[1] = 0.6521451549;
        fPoints(2,0) = 0.3399810436;        fWeights[2] = 0.6521451549;
        fPoints(3,0) = 0.8611363116;        fWeights[3] = 0.3478548451;
    }
    if (fOrder == 8 || fOrder == 9) {            
        fPoints(0,0) = -0.9061798459;        fWeights[0] = 0.2369268851;
        fPoints(1,0) = -0.5384693101;        fWeights[1] = 0.4786286705;
        fPoints(2,0) = 0.0000000000;        fWeights[2] = 0.5688888889;
        fPoints(3,0) = 0.5384693101;        fWeights[3] = 0.4786286705;
        fPoints(4,0) = 0.9061798459;        fWeights[4] = 0.2369268851;
    }
    if (fOrder == 10 || fOrder == 11) {        
        fPoints(0,0) = -0.9324695142;        fWeights[0] = 0.1713244924;
        fPoints(1,0) = -0.6612093865;        fWeights[1] = 0.360761573;
        fPoints(2,0) = -0.2386191861;        fWeights[2] = 0.4679139346;
        fPoints(3,0) = 0.2386191861;        fWeights[3] = 0.4679139346;
        fPoints(4,0) = 0.6612093865;        fWeights[4] = 0.360761573;
        fPoints(5,0) = 0.9324695142;        fWeights[5] = 0.1713244924;
    }
    if (fOrder == 12 || fOrder == 13) {        
        fPoints(0,0) = -0.9491079123;        fWeights[0] = 0.1294849662;
        fPoints(1,0) = -0.7415311856;        fWeights[1] = 0.2797053915;
        fPoints(2,0) = -0.4058451514;        fWeights[2] = 0.3818300505;
        fPoints(3,0) = 0.00000000000;        fWeights[3] = 0.4179591837;
        fPoints(4,0) = 0.4058451514;        fWeights[4] = 0.3818300505;
        fPoints(5,0) = 0.7415311856;        fWeights[5] = 0.2797053915;
        fPoints(6,0) = 0.9491079123;        fWeights[6] = 0.1294849662;
    }
    if (fOrder == 14 || fOrder == 15) {
        fPoints(0,0) = -0.9602898565;        fWeights[0] = 0.1012285363;
        fPoints(1,0) = -0.7966664774;        fWeights[1] = 0.2223810345;
        fPoints(2,0) = -0.5255324099;        fWeights[2] = 0.3137066459;
        fPoints(3,0) = -0.1834346425;        fWeights[3] = 0.3626837834;
        fPoints(4,0) = 0.1834346425;        fWeights[4] = 0.3626837834;
        fPoints(5,0) = 0.5255324099;        fWeights[5] = 0.3137066459;
        fPoints(6,0) = 0.7966664774;        fWeights[6] = 0.2223810345;
        fPoints(7,0) = 0.9602898565;        fWeights[7] = 0.1012285363;
    }
    if (fOrder == 16 || fOrder == 17) {            
        fPoints(0,0) = -0.9681602395;        fWeights[0] = 0.08127438836;
        fPoints(1,0) = -0.8360311073;        fWeights[1] = 0.1806481607;
        fPoints(2,0) = -0.6133714327;        fWeights[2] = 0.2606106964;
        fPoints(3,0) = -0.3242534234;        fWeights[3] = 0.312347077;
        fPoints(4,0) = 0.0000000000;        fWeights[4] = 0.330239355;
        fPoints(5,0) = 0.3242534234;        fWeights[5] = 0.312347077;
        fPoints(6,0) = 0.6133714327;        fWeights[6] = 0.2606106964;
        fPoints(7,0) = 0.8360311073;        fWeights[7] = 0.1806481607;
        fPoints(8,0) = 0.9681602395;        fWeights[8] = 0.08127438836;
    }
    if (fOrder == 18 || fOrder == 19) {            
        fPoints(0,0) = -0.9739065285;        fWeights[0] = 0.06667134431;
        fPoints(1,0) = -0.8650633667;        fWeights[1] = 0.1494513492;
        fPoints(2,0) = -0.6794095683;        fWeights[2] = 0.2190863625;
        fPoints(3,0) = -0.4333953941;        fWeights[3] = 0.2692667193;
        fPoints(4,0) = -0.148874339;        fWeights[4] = 0.2955242247;
        fPoints(5,0) = 0.148874339;        fWeights[5] = 0.2955242247;
        fPoints(6,0) = 0.4333953941;        fWeights[6] = 0.2692667193;
        fPoints(7,0) = 0.6794095683;        fWeights[7] = 0.2190863625;
        fPoints(8,0) = 0.8650633667;        fWeights[8] = 0.1494513492;
        fPoints(9,0) = 0.9739065285;        fWeights[9] = 0.06667134431;
    }
    if (fOrder > 19) {
        VecDouble co(npoints);
        VecDouble w(npoints);

        this->gauleg(-1, 1, co, w);

        for (int i = 0; i < npoints; i++) {
            fPoints(i, 0) = co[i];
            fWeights[i] = w[i];
        }
        co.clear();
        w.clear();
    }
}