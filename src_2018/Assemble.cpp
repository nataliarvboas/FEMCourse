/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Assemble.h"
#include "GeoElement.h"

Assemble::Assemble() {
}

Assemble::Assemble(CompMesh *mesh) {
    cmesh = mesh;
}

Assemble::Assemble(const Assemble &copy) {
    cmesh = copy.cmesh;
}

Assemble &Assemble::operator=(const Assemble &copy) {
    cmesh = copy.cmesh;
    return *this;
}

void Assemble::SetMesh(CompMesh *mesh) {
    cmesh = mesh;
}

int64_t Assemble::NEquations() {
    int64_t neq = 0;
    int64_t i, ncon = cmesh->GetDOFVec().size();
    for (i = 0; i < ncon; i++) {
        DOF dof = cmesh->GetDOF(i);
        int dofsize = dof.GetNShape() * dof.GetNState();
        neq += dofsize;
    }
    return neq;
}

void Assemble::OptimizeBandwidth() {
}

void Assemble::Compute(Matrix &globmat, Matrix &rhs) {
    int64_t nelem = cmesh->GetGeoMesh()->NumElements();

    globmat.Zero();
    rhs.Zero();

    for (int64_t el = 0; el < nelem; el++) {
        CompElement *cel = cmesh->GetElement(el);
        GeoElement *gel = cmesh->GetGeoMesh()->Element(el);

        int64_t nnodes_el = cel->GetGeoElement()->NNodes();

        Matrix ek(nnodes_el, nnodes_el);
        Matrix ef(nnodes_el, 1);

        ek.Zero();
        ef.Zero();

        cel->CalcStiff(ek, ef);

        for (int64_t i = 0; i < nnodes_el; i++) {
            int64_t IG = gel->NodeIndex(i);
            rhs(IG, 0) += ef(i, 0);
            
            for (int64_t j = 0; j < nnodes_el; j++) {
                int64_t JG = gel->NodeIndex(j);
                globmat(IG, JG) += ek(i, j);
            }
        }
    }
}