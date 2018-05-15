/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include <vector> 
#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTriangle.h"
#include "IntRuleTetrahedron.h"
#include <cmath>
#include <iomanip>

using namespace std;

IntRule::IntRule() {

}

IntRule::IntRule(int order) {
    fOrder = order;
}

IntRule::~IntRule() {

}

IntRule::IntRule(const IntRule& copy) {
    fOrder = copy.fOrder;
}

IntRule &IntRule::operator=(const IntRule &cp) {
    fOrder = cp.fOrder;
    return *this;
}

int IntRule::NPoints() const {
    return fWeights.size();
}

void IntRule::Print(std::ostream &out) const {
    VecDouble co(fPoints.Cols());
    double w;

    for (int i = 0; i < NPoints(); i++) {
        Point(i, co, w);

        int dim = co.size();

        out << "ponto: " << i << endl;

        for (int j = 0; j < dim; j++) {
            out << "coord: " << setprecision(10) << co[j] << endl;
        }
        out << "peso: " << setprecision(10) << w << endl << endl;
    }
}

void IntRule::Point(int p, VecDouble& co, double& w) const {

    int dim = co.size();

    for (int i = 0; i < dim; i++) {
        co[i] = fPoints.GetVal(p, i);
    }
    w = fWeights[p];
}

