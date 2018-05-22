/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
#include "CompElement.h"

Poisson::Poisson() {
}

Poisson::Poisson(Matrix &perm) {
    permeability = perm;
}

Poisson::Poisson(const Poisson &copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
}

Poisson &Poisson::operator=(const Poisson &copy) {
    permeability = copy.permeability;
    forceFunction = copy.forceFunction;
    return *this;
}

Poisson *Poisson::Clone() const {
    return new Poisson(*this);
}

Poisson::~Poisson() {
}

Matrix Poisson::GetPermeability() const {
    return permeability;
}

void Poisson::SetPermeability(const Matrix &perm) {
    permeability = perm;
}

int Poisson::NState() const {
    return 1;
}

void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    VecDouble phi = data.phi;
    Matrix dphi = data.dphidx;
    Matrix axes = data.axes;

    int nshape = phi.size();
    int dim = dphi.Rows();

    std::vector<double> du(dim);
    std::vector<double> dv(dim);


    Matrix perm(dim, dim);
    std::vector<double> result(nshape);

    perm = this->GetPermeability();
    //force = this->GetForceFunction();

    for (int i = 0; i < nshape; i++) {
        dv[0] = dphi(0, i) * axes(0, 0) + dphi(1, i) * axes(1, 0);
        dv[1] = dphi(0, i) * axes(0, 1) + dphi(1, i) * axes(1, 1);

        for (int j = 0; j < nshape; j++) {
            du[0] = dphi(0, j) * axes(0, 0) + dphi(1, j) * axes(1, 0);
            du[1] = dphi(0, j) * axes(0, 1) + dphi(1, j) * axes(1, 1);

            const int posI = 2 * i;
            const int posJ = 2 * j;

            EK(posI, posJ) += du[0] * dv[0] * perm(0, 0) * weight + du[0] * dv[1] * perm(1, 0) * weight;
            EK(posI, posJ + 1) += du[0] * dv[0] * perm(0, 1) * weight + du[0] * dv[1] * perm(1, 1) * weight;
            EK(posI + 1, posJ) += du[1] * dv[0] * perm(0, 0) * weight + du[1] * dv[1] * perm(1, 0) * weight;
            EK(posI + 1, posJ + 1) += du[1] * dv[0] * perm(0, 1) * weight + du[1] * dv[1] * perm(1, 1) * weight;
        }
        for (int k = 0; k < dim; k++) {
            //EF(i, 0) += phi[i] * result[i] * weight;
            EF(i, 0) += phi[i] * weight;
        }
    }
}

std::vector<double> Poisson::PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const {
}
