/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "CompMesh.h"

CompMesh::CompMesh() {
}

CompMesh::CompMesh(const CompMesh &copy) {
    compelements = copy.compelements;
    dofs = copy.dofs;
    mathstatements = copy.mathstatements;
}

CompMesh::~CompMesh() {
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

void CompMesh::Resequence() {
}

void CompMesh::Resequence(VecInt &DOFindices) {
}
    