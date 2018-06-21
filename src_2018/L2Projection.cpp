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

L2Projection::L2Projection(int bctype, int materialid, Matrix &proj, Matrix Val1, Matrix Val2) {
    projection = proj;
    BCType = bctype;
    BCVal1 = Val1;
    BCVal2 = Val2;
    this->SetMatID(materialid);
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;

}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;
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
    return 3;
}

void L2Projection::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    int nstate = this->NState();
    int nshape = data.phi.size();

    VecDouble result(data.x.size());
    Matrix deriv(data.x.size(), data.x.size());

    SolutionExact(data.x, result, deriv);

    switch (this->GetBCType()) {

        case 0:
        {
            for (int iv = 0; iv < nstate; iv++) {
                for (int in = 0; in < nshape; in++) {
                    EF(nstate * in + iv, 0) += MathStatement::gBigNumber * result[iv] * data.phi[in] * weight;
                    for (int jn = 0; jn < nshape; jn++) {
                        EK(nstate * in + iv, nstate * jn + iv) += MathStatement::gBigNumber * data.phi[in] * data.phi[jn] * weight;
                    }
                }
            }
            break;
        }

        case 1:
        {
            for (int iv = 0; iv < nstate; iv++) {
                for (int in = 0; in < nshape; in++) {
                    EF(nstate * in + iv, 0) += Val2()(iv, 0) * data.phi[in] * weight;
                }
            }
            break;
        }

        default:
        {
            std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
        }
    }

    //    for (int i = 0; i < nshape; i++) {
    //        for (int ivi = 0; ivi < nstate; ivi++) {
    //            const int posI = nstate * i + ivi;
    //            EF(posI, 0) += weight * data.phi[i] * result[ivi] * MathStatement::gBigNumber;
    //        }
    //        for (int j = 0; j < nshape; j++) {
    //            for (int ivi = 0; ivi < nstate; ivi++) {
    //                const int posI = nstate * i + ivi;
    //                const int posJ = nstate * j + ivi;
    //                EK(posI, posJ) += weight * data.phi[i] * data.phi[j] * MathStatement::gBigNumber;
    //            }
    //        }
    //    }
}

int L2Projection::NEvalErrors() const {
    return 3;
}

void L2Projection::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    return;
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
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->NState();

}

std::vector<double> L2Projection::PostProcessSolution(const IntPointData &integrationpointdata, const int var) const {
    //    DebugStop();
    VecDouble vec(2, 0);
    return vec;

}