/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
#include "PostProcess.h"
//#include "CompElement.h"
#include <functional>
#include <string.h>

Poisson::Poisson() {
}

Poisson::Poisson(int materialid, Matrix &perm) {
    permeability = perm;
    this->SetMatID(materialid);
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

int Poisson::NEvalErrors() const {
    return 3;
}

int Poisson::NState() const {
    return 2;
}

int Poisson::VariableIndex(const PostProcVar var) const {
    if (var == ESol) return ESol;
    if (var == EDSol) return EDSol;
    if (var == EFlux) return EFlux;
    if (var == EForce) return EForce;
    if (var == ESolExact) return ESolExact;
    if (var == EDSolExact) return EDSolExact;
}

Poisson::PostProcVar Poisson::VariableIndex(const std::string &name) {
    if (!strcmp("Sol", name.c_str())) return ESol;
    if (!strcmp("DSol", name.c_str())) return EDSol;
    if (!strcmp("Flux", name.c_str())) return EFlux;
    if (!strcmp("Force", name.c_str())) return EForce;
    if (!strcmp("SolExact", name.c_str())) return ESolExact;
    if (!strcmp("DSolExact", name.c_str())) return EDSolExact;
}

int Poisson::NSolutionVariables(const PostProcVar var) {


}

void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
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

void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    VecDouble phi = data.phi;
    Matrix dphi = data.dphidx;
    Matrix axes = data.axes;

    int nshape = phi.size();
    int dim = this->Dimension();
    int nstate = this->NState();

    Matrix perm(dim, dim);
    std::function<void(const VecDouble &co, VecDouble & result) > force;

    perm = this->GetPermeability();
    force = this->GetForceFunction();

    VecDouble res(data.x.size());
    force(data.x, res);

    Matrix nvec(nstate, nstate, 0);
    for (int i = 0; i < nstate; i++) {
        nvec(i, i) = 1;
    }

    int id = 0;
    int inormal = 0;
    VecDouble normalid(nshape*nstate, 0.);
    VecDouble shapeid(nshape*nstate, 0.);
    for (int i = 0; i < nshape; i++) {
        inormal = 0;
        for (int s = 0; s < nstate; s++) {
            normalid[id] = inormal;
            shapeid[id] = i;
            id++;
            inormal++;
        }
    }

    Matrix dphi_nvec(nstate, nstate);
    std::vector<Matrix> grad(nshape * nstate);
    int i, j, k, nvecid;
    for (i = 0; i < nshape * nstate; i++) {
        nvecid = normalid[i];
        for (j = 0; j < nstate; j++) {
            for (k = 0; k < nstate; k++) {
                dphi_nvec(j, k) = dphi(k, shapeid[i]) * nvec(j, nvecid);
            }
        }
        grad[i] = dphi_nvec;
    }

    Matrix mult(nstate, nstate);
    double inner;

    for (int i = 0; i < nshape; i++) {
        for (int k = 0; k < nstate; k++) {
            int posI = nstate * i + k;
            EF(posI, 0) += phi[i] * res[k] * weight;
        }
    }

    for (int i = 0; i < nshape * nstate; i++) {
        for (int j = 0; j < nshape * nstate; j++) {
            mult = grad[j] * perm;
            inner = this->Inner(grad[i], mult);
            EK(i, j) += inner * weight;
        }
    }

    //        for (int i = 0; i < nshape; i++) {
    //            const int posI = 2 * i;
    //    
    //            EF(posI, 0) += phi[i] * res[0] * weight;
    //            EF(posI + 1, 0) += phi[i] * res[1] * weight;
    //    
    //            dv[0] = dphi(0, i) * axes(0, 0) + dphi(1, i) * axes(1, 0);
    //            dv[1] = dphi(0, i) * axes(0, 1) + dphi(1, i) * axes(1, 1);
    //    
    //            for (int j = 0; j < nshape; j++) {
    //                const int posJ = 2 * j;
    //                du[0] = dphi(0, j) * axes(0, 0) + dphi(1, j) * axes(1, 0);
    //                du[1] = dphi(0, j) * axes(0, 1) + dphi(1, j) * axes(1, 1);
    //    
    //    
    //                EK(posI, posJ) += du[0] * dv[0] * perm(0, 0) * weight + du[0] * dv[1] * perm(1, 0) * weight;
    //                EK(posI, posJ + 1) += du[0] * dv[0] * perm(0, 1) * weight + du[0] * dv[1] * perm(1, 1) * weight;
    //                EK(posI + 1, posJ) += du[1] * dv[0] * perm(0, 0) * weight + du[1] * dv[1] * perm(1, 0) * weight;
    //                EK(posI + 1, posJ + 1) += du[1] * dv[0] * perm(0, 1) * weight + du[1] * dv[1] * perm(1, 1) * weight;
    //            }
    //        }
}

std::vector<double> Poisson::PostProcessSolution(const IntPointData &integrationpointdata, const int var) const {

}

double Poisson::Inner(Matrix &S, Matrix & T) const {
    double inner = 0;
    for (int i = 0; i < S.Rows(); i++) {
        for (int j = 0; j < S.Cols(); j++) {
            inner += S(i, j) * T(i, j);
        }
    }
    return inner;
}
