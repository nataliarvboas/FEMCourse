/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"

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

void Poisson::Contribute(IntPointData &integrationpointdata, double weight , Matrix &EK, Matrix &EF) const{
}

void Poisson::PostProcess(IntPointData &integrationpointdata, const std::string &variable, VecDouble &postprocvalue) const {
}
