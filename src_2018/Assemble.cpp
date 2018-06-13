/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Assemble.h"
//#include "GeoElement.h"
#include "CompMesh.h"
#include "GeoMesh.h"
//#include <algorithm> 

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
    int64_t ne = this->NEquations();

    globmat.Resize(ne, ne);
    rhs.Resize(ne, 1);

    globmat.Zero();
    rhs.Zero();

    for (int64_t el = 0; el < nelem; el++) {
        CompElement *cel = cmesh->GetElement(el);

        Matrix ek, ef;
        ek.Zero();
        ef.Zero();
        cel->CalcStiff(ek, ef);
        
        ek.Print();
        ef.Print();

        int ndof = cel->NDOF();

        VecInt iglob(0);
        int ni = 0;
        for (int i = 0; i < ndof; i++) {
            int dofindex = cel->GetDOFIndex(i);
            DOF dof = cmesh->GetDOF(dofindex);
            for (int j = 0; j < dof.GetNShape() * dof.GetNState(); j++) {
                iglob.resize(ni + 1);
                iglob[ni] = dof.GetFirstEquation() + j;
                ni++;
            }
        }

        for (int64_t i = 0; i < ek.Rows(); i++) {
            int64_t IG = iglob[i];
            rhs(IG, 0) += ef(i, 0);

            for (int64_t j = 0; j < ek.Rows(); j++) {
                int64_t JG = iglob[j];
                globmat(IG, JG) += ek(i, j);
            }
        }
    }
}
