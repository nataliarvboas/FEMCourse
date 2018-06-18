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
    if (var == ESol) return this->Dimension();
    if (var == EDSol) return this->Dimension();
    if (var == EFlux) return this->Dimension();
    if (var == EForce) return this->Dimension();
    if (var == ESolExact) return this->Dimension();
    if (var == EDSolExact) return this->Dimension();
    else {
        std::cout << "variable not implemented" << std::endl;
    }

}

void Poisson::ContributeError(IntPointData &data, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const {
    errors.resize(NEvalErrors());

    for (int i = 0; i < this->NEvalErrors(); i++) {
        errors[i] = 0;
    }
    Matrix dudx = data.dsoldx;
    Matrix gradu(3, 1);
    Matrix axes = data.axes;
    this->Axes2XYZ(dudx, gradu, axes);

    VecDouble sol(1), dsol(3, 0.);
    VecDouble u = data.solution;
    double diff;

    for (int i = 0; i < this->NState(); i++) {
        diff += (u[i] - u_exact[i]);
    }


    errors[1] = diff*diff;

    errors[2] = 0.;

    for (int i = 0; i < this->Dimension(); i++) {
        for (int j = 0; j < this->NState(); j++) {
            diff = (gradu(i, j) - du_exact(i, j));
            errors[2] += diff*diff;
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

    Matrix dphi_(dim, 1);
    std::vector<Matrix> grad(nshape);
    int i, j, k;
    for (i = 0; i < nshape; i++) {
        for (j = 0; j < dim; j++) {
            dphi_(j, 0) = dphi(j, i);
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
                int posI = nstate * i + ivi;
                int posJ = nstate * j + ivi;
                Matrix gradi_T;
                grad[i].Transpose(gradi_T);
                EK(posI, posJ) += (gradi_T * grad[j])(0, 0) * weight;
            }
        }
    }
}

std::vector<double> Poisson::PostProcessSolution(const IntPointData &data, const int var) const {
    VecDouble sol = data.solution;
    int rows = data.dsoldx.Rows();
    int cols = data.dsoldx.Cols();
    Matrix gradu(rows, cols);
    gradu = data.dsoldx;

    int dim = this->Dimension();
    VecDouble Solout(0);

    switch (var) {
        case 0: //None
        {
            break;
        }

        case 1: //ESol
        {
            Solout.resize(dim);
            Solout[0] = sol[0];
            Solout[1] = sol[1];
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
            Solout.resize(dim);
            VecDouble result(dim);
            this->forceFunction(data.x, result);
            Solout[0] = result[0];
            Solout[1] = result[1];
        }
            break;

        case 5: //ESolExact
        {
            Solout.resize(dim);
            VecDouble sol(dim);
            Matrix dsol(dim, dim);
            this->SolutionExact(data.x, sol, dsol);

            Solout[0] = sol[0];
            Solout[1] = sol[1];

        }
            break;
        case 6: //EDSolExact
        {
            Solout.resize(rows * cols);
            VecDouble sol(dim);
            Matrix dsol(dim, dim);
            this->SolutionExact(data.x, sol, dsol);

            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    Solout[i * cols + j] = dsol(i, j);
                }
            }
        }
            break;


        default:
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return Solout;

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
