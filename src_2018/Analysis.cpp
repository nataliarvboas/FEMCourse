/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Analysis.h"
#include "Assemble.h"
#include "CompMesh.h"

Analysis::Analysis() {
    cmesh = 0;
    Solution = 0;
    GlobalSystem = 0;
    RightHandSide = 0;

}

Analysis::Analysis(const Analysis &cp) {
    cmesh = cp.cmesh;
    Solution = cp.Solution;
    GlobalSystem = cp.GlobalSystem;
    RightHandSide = cp.RightHandSide;
}

Analysis &Analysis::operator=(const Analysis &cp) {
    cmesh = cp.cmesh;
    Solution = cp.Solution;
    GlobalSystem = cp.GlobalSystem;
    RightHandSide = cp.RightHandSide;
    return *this;
}

Analysis::Analysis(CompMesh *mesh) {
    cmesh = mesh;
}

void Analysis::SetMesh(CompMesh *mesh) {
    cmesh = mesh;
}

CompMesh *Analysis::Mesh() const {
    return cmesh;
}

void Analysis::RunSimulation() {

    double nnodes = cmesh->GetElementVec().size();
    Matrix this->GlobalSystem.Resize(nnodes, nnodes, 0);
    GlobalSystem.Zero();
    Matrix RightHandSide(nnodes, 1, 0);

    Assemble assemb = new Assemble;

    assemb->Compute(GlobalSystem, RightHandSide);

    GlobalSystem.Print();
    RightHandSide.Print();
}

//void Analysis::PostProcess(std::string &filename, PostProcess &defPostProc) const {
//
//}
