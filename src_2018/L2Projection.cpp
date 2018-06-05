/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"

L2Projection::L2Projection() {
}

L2Projection::L2Projection(int materialid, Matrix &perm) {
    projection = perm;
    this->SetMatID(materialid);
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    return *this;
}

L2Projection *L2Projection::Clone() const {
    return new L2Projection(*this);
}

L2Projection::~L2Projection() {
}

Matrix L2Projection::GetProjectionMatrix() const {
    return projection;
}

void L2Projection::SetProjectionMatrix(const Matrix &proj) {
    projection = proj;
}

int L2Projection::NState() const {
    return 1;

}

void L2Projection::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    int nvars = this->NState();
    int nshape = data.phi.size();
    for (int i = 0; i < nshape; i++) {
        for (int j = 0; j < nshape; j++) {
            for (int ivi = 0; ivi < nvars; ivi++) {
                const int posI = nvars * i + ivi;
                const int posJ = nvars * j + ivi;
                EK(posI, posJ) += weight * data.phi[i] * data.phi[j];
            }
        }
        for (int ivi = 0; ivi < nvars; ivi++) {
            const int posI = nvars * i + ivi;
            EF(posI, 0) += weight * data.phi[i] * data.solution[ivi];
        }
    }
}

int L2Projection::NEvalErrors() const {
}

void L2Projection::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
}

int L2Projection::VariableIndex(const std::string &name){
}

int L2Projection::NSolutionVariables(const PostProcVar var){
}

std::vector<double> L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const PostProcVar var) const {

}