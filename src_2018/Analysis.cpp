/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Analysis.h"
#include "Assemble.h"

Analysis::Analysis() {
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
    Assemble *assemb = new Assemble();

    int64_t neq = assemb->NEquations();

    if (neq > 20000) {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }

    Assemble();

    if (neq > 20000) {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
}

void Analysis::PostProcess(std::string &filename, PostProcess &defPostProc) const {
    
}
