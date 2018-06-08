/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"
#include "PostProcess.h"

L2Projection::L2Projection() {
}

L2Projection::L2Projection(int materialid, Matrix &perm) {
    projection = perm;
    this->SetMatID(materialid);
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
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
    return 2;
}

void L2Projection::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    int nvars = this->NState();
    int nshape = data.phi.size();

    VecDouble result(data.x.size());
    Matrix deriv(data.x.size(), data.x.size());

    SolutionExact(data.x, result, deriv);

    for (int i = 0; i < nshape; i++) {
        const int posI = 2 * i;
        EF(posI, 0) += weight * data.phi[i] * result[0] * MathStatement::gBigNumber;
        EF(posI + 1, 0) += weight * data.phi[i] * result[1] * MathStatement::gBigNumber;
        for (int j = 0; j < nshape; j++) {
            const int posJ = 2 * j;
            EK(posI, posJ) += weight * data.phi[i] * data.phi[j] * MathStatement::gBigNumber;
            EK(posI+ 1, posJ + 1) += weight * data.phi[i] * data.phi[j] * MathStatement::gBigNumber;
        }
    }
}

int L2Projection::NEvalErrors() const {
}

void L2Projection::ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
}

int L2Projection::VariableIndex(const PostProcVar var) const {
}
//

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string &name) {
}

int L2Projection::NSolutionVariables(const PostProcVar var) {
}

std::vector<double> L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var) const {

}