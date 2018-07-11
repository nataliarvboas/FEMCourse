/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "Poisson.h"
#include "PostProcess.h"
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

int Poisson::VariableIndex(const PostProcVar var) const {
    if (var == ENone) return ENone;
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
    else {
        std::cout << "variable not implemented" << std::endl;
    }
}

int Poisson::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->Dimension();
    if (var == EFlux) return this->Dimension();
    if (var == EForce) return this->NState();
    if (var == ESolExact) return this->NState();
    if (var == EDSolExact) return this->Dimension();
    else {
        std::cout << "variable not implemented" << std::endl;
    }

}

void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    errors.resize(NEvalErrors(), 0);
    Matrix gradu;
    Matrix axes = data.axes;

    VecDouble u = data.solution;
    Matrix dudx = data.dsoldx;

    this->Axes2XYZ(dudx, gradu, axes);

    double diff = 0.0;
    for (int i = 0; i < this->NState(); i++) {
        diff = (u[i] - u_exact[i]);
        errors[0] += diff*diff;
    }

    errors[1] = 0.;
    int dim = this->Dimension();
    int nstate = this->NState();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < nstate; j++) {
            diff = (gradu(i, j) - du_exact(i, j));
            errors[1] += diff*diff;
        }

    }
    errors[2] = errors[0] + errors[1];
    u.clear();
}

void Poisson::Contribute(IntPointData &data, double weight, Matrix &EK, Matrix &EF) const {
    VecDouble phi = data.phi;
    Matrix dphi = data.dphidx;
    Matrix axes = data.axes;
    Matrix dphi2;

    this->Axes2XYZ(dphi, dphi2, axes);

    int nshape = phi.size();
    int nstate = this->NState();
    int dim = dphi.Rows();

    Matrix perm(dim, dim);
    std::function<void(const VecDouble &co, VecDouble & result) > force;

    perm = this->GetPermeability();
    force = this->GetForceFunction();

    VecDouble res(data.x.size());
    force(data.x, res);

    Matrix dphi_(dim, 1);
    std::vector<Matrix> grad(nshape);
    int i, j, k;
    for (i = 0; i < nshape; i++) {
        for (j = 0; j < dim; j++) {
            dphi_(j, 0) = dphi2(j, i);
        }
        grad[i] = dphi_;
    }

    int ivi = 0;
    for (i = 0; i < nshape; i++) {
        for (ivi = 0; ivi < nstate; ivi++) {
            int posI = nstate * i + ivi;
            EF(posI, 0) += phi[i] * res[ivi] * weight;
        }
        for (j = 0; j < nshape; j++) {
            for (ivi = 0; ivi < nstate; ivi++) {
                for (int ivj = 0; ivj < nstate; ivj++) {
                    int posI = nstate * i + ivi;
                    int posJ = nstate * j + ivj;
                    Matrix gradi_T;
                    grad[i].Transpose(gradi_T);
                    EK(posI, posJ) += (gradi_T * grad[j])(0, 0) * perm(ivi, ivj) * weight;
                }
            }
        }
    }
    phi.clear();
    res.clear();
    grad.clear();
}

void Poisson::PostProcessSolution(const IntPointData &data, const int var, VecDouble &Solout) const {
    VecDouble sol = data.solution;
    int solsize = sol.size();
    int rows = data.dsoldx.Rows();
    int cols = data.dsoldx.Cols();
    Matrix gradu(rows, cols);
    gradu = data.dsoldx;

    int nstate = this->NState();

    switch (var) {
        case 0: //None
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }

        case 1: //ESol
        {
            Solout.resize(nstate);
            for (int i = 0; i < nstate; i++) {
                Solout[i] = sol[i];
            }
        }
            break;

        case 2: //EDSol
        {
            Solout.resize(rows * cols);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    Solout[i * cols + j] = gradu(i, j);
                }
            }

        }
            break;
        case 3: //EFlux
        {

        }
            break;

        case 4: //EForce
        {
            Solout.resize(nstate);
            VecDouble result(nstate);
            this->forceFunction(data.x, result);
            for (int i = 0; i < nstate; i++) {
                Solout[i] = result[i];
            }
            result.clear();
        }
            break;

        case 5: //ESolExact
        {
            Solout.resize(nstate);
            VecDouble sol(nstate);
            Matrix dsol(nstate, nstate);
            this->SolutionExact(data.x, sol, dsol);
            for (int i = 0; i < nstate; i++) {
                Solout[i] = sol[i];
            }
            sol.clear();

        }
            break;
        case 6: //EDSolExact
        {
            Solout.resize(rows * cols);
            VecDouble sol(nstate);
            Matrix dsol(nstate, nstate);
            this->SolutionExact(data.x, sol, dsol);

            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    Solout[i * cols + j] = dsol(i, j);
                }
            }
            sol.clear();
        }
            break;


        default:
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    sol.clear();
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