/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"
#include "PostProcess.h"
#include "tpanic.h"
#include <string.h>

L2Projection::L2Projection() {
}

L2Projection::L2Projection(int bctype, int materialid, Matrix &perm) {
    projection = perm;
    BCType = bctype;
    this->SetMatID(materialid);
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;

}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
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
        for (int ivi = 0; ivi < nvars; ivi++) {
            const int posI = nvars * i + ivi;
            EF(posI, 0) += weight * data.phi[i] * result[ivi] * MathStatement::gBigNumber;
        }
        for (int j = 0; j < nshape; j++) {
            for (int ivi = 0; ivi < nvars; ivi++) {
                const int posI = nvars * i + ivi;
                const int posJ = nvars * j + ivi;
                EK(posI, posJ) += weight * data.phi[i] * data.phi[j] * MathStatement::gBigNumber;
            }
        }
    }
}

int L2Projection::NEvalErrors() const {
    return 3;
}

void L2Projection::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    errors.resize(NEvalErrors());

    for (int i = 0; i < this->NEvalErrors(); i++) {
        errors[i] = 0;
    }
    Matrix dudx = data.dsoldx;
    Matrix gradu(3, 1);
    Matrix axes = data.axes;
    this->Axes2XYZ(dudx, gradu, axes, true);

    VecDouble sol(1), dsol(3, 0.);
    VecDouble u = data.solution;
    double diff;
    
    for (int i = 0; i < u.size(); i++) {
        diff = (u[0] - u_exact[0]);
    }


    errors[1] = diff*diff;

    errors[2] = 0.;

    for (int i = 0; i < gradu.Rows(); i++) {
        for (int j = 0; j < gradu.Cols(); j++) {
            diff += (gradu(i, j) - du_exact(i, j));
        }
    }
    errors[0] = errors[1] + errors[2];
}

int L2Projection::VariableIndex(const PostProcVar var) const {
    if (var == ESol) return ESol;
    if (var == EDSol) return EDSol;
}

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string & name) {
    if (!strcmp("Solution", name.c_str())) return ESol;
    if (!strcmp("Derivative", name.c_str())) return EDSol;
}

int L2Projection::NSolutionVariables(const PostProcVar var) {
    //    if (var == ESol) return this->NState();
    //    if (var == EDSol) return this->NState(); //??

}

std::vector<double> L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var) const {
    DebugStop();
}