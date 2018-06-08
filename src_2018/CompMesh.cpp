/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "CompMesh.h"
//#include "MathStatement.h"
#include "GeoMesh.h"

CompMesh::CompMesh() : geomesh(0), compelements(0), dofs(0), mathstatements(0), solution(0) {

}

CompMesh::CompMesh(const CompMesh &copy) {
    compelements = copy.compelements;
    dofs = copy.dofs;
    mathstatements = copy.mathstatements;
}

CompMesh::CompMesh(GeoMesh *gmesh) {
    this->SetGeoMesh(gmesh);
    this->SetNumberElement(gmesh->NumElements());
}

CompMesh::~CompMesh() {
}

GeoMesh *CompMesh::GetGeoMesh() const {
    return geomesh;
}

void CompMesh::SetGeoMesh(GeoMesh *gmesh) {
    geomesh = gmesh;
}

void CompMesh::SetNumberElement(int64_t nelem) {
    compelements.resize(nelem);
}

void CompMesh::SetNumberDOF(int64_t ndof) {
    dofs.resize(ndof);
}

void CompMesh::SetNumberMath(int nmath) {
    mathstatements.resize(nmath);
}

void CompMesh::SetElement(int64_t elindex, CompElement *cel) {
    compelements[elindex] = cel;
}

void CompMesh::SetDOF(int64_t index, const DOF &dof) {
    dofs[index] = dof;
}

void CompMesh::SetMathStatement(int index, MathStatement *math) {
    mathstatements[index] = math;
}

DOF &CompMesh::GetDOF(int64_t dofindex) {
    return dofs[dofindex];
}

CompElement *CompMesh::GetElement(int64_t elindex) const {
    return compelements[elindex];
}

MathStatement *CompMesh::GetMath(int matindex) const {
    return mathstatements[matindex];
}

std::vector<CompElement *> CompMesh::GetElementVec() const {
    return compelements;
}

std::vector<DOF> CompMesh::GetDOFVec() const {
    return dofs;
}

std::vector<MathStatement *> CompMesh::GetMathVec() const {
    return mathstatements;
}

void CompMesh::SetElementVec(const std::vector<CompElement *> &vec) {
    int c = vec.size();

    for (int i = 0; i < c; i++) {
        compelements[i] = vec[i];
    }
}

void CompMesh::SetDOFVec(const std::vector<DOF> &dofvec) {
    int c = dofvec.size();

    for (int i = 0; i < c; i++) {
        dofs[i] = dofvec[i];
    }
}

void CompMesh::SetMathVec(const std::vector<MathStatement *> &mathvec) {
    int c = mathvec.size();

    for (int i = 0; i < c; i++) {
        mathstatements[i] = mathvec[i];
    }
}

void CompMesh::AutoBuild() {
    int64_t i;
    int64_t nelem = this->GetElementVec().size();

    for (i = 0; i < nelem; i++) {
        GeoElement *gel = this->GetGeoMesh()->Element(i);
        CompElement *cel = gel->CreateCompEl(this, i);
        this->SetElement(i, cel);
    }
    this->Resequence();
}

void CompMesh::Resequence() {
    int64_t ndof = this->GetNumberDOF();
    int64_t fe = 0;
    this->GetDOF(0).SetFirstEquation(0);

    for (int i = 1; i < ndof; i++) {
        int nshape = this->GetDOF(i).GetNShape();
        int nstate = this->GetDOF(i).GetNState();
        int result = nshape * nstate;
        fe += result;
        this->GetDOF(i).SetFirstEquation(fe);
    }
}

void CompMesh::Resequence(VecInt &DOFindices) {
}

std::vector<double> &CompMesh::Solution() {
    return solution;
}

void CompMesh::LoadSolution(std::vector<double> &Sol) {
    solution.resize(Sol.size());
    for (int64_t i = 0; i < Sol.size(); i++) {
        solution[i] = Sol[i];
    }
}
