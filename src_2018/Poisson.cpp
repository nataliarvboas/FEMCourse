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
    //return new Poisson(*this);
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
    int nshape = data.phi.size();
    CompElement *comp;
    int dim = comp->Dimension();    
    Matrix perm= this->GetPermeability();
    
    for (int in = 0; in < nshape; in++) {
        int d;
        for (d = 0; d < dim; d++) {
            for (int jn = 0; jn < nshape; jn++) {
                EK(in, jn) += weight * data.dphidksi(d, in) * data.dphidksi(d, jn) * perm(in, jn);
            }
            EF(in, 0) += -weight * data.phi[in] * data.solution[in];
        }
    }
}

std::vector<double> Poisson::PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const {
}
