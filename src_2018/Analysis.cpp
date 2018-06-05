/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Analysis.h"
#include "Assemble.h"
#include "CompMesh.h"
#include "GeoMesh.h"

Analysis::Analysis() {
    cmesh = 0;
    Solution.Resize(0, 0);
    GlobalSystem.Resize(0, 0);
    RightHandSide.Resize(0, 0);

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

Analysis::~Analysis() {
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
    double nnodes = cmesh->GetGeoMesh()->NumNodes();
    GlobalSystem.Resize(nnodes, nnodes);
    RightHandSide.Resize(nnodes, 1);
    Matrix F(nnodes, 1);

    Assemble assemb(cmesh);
    assemb.Compute(GlobalSystem, RightHandSide);

    std::cout << "\nGlobal Stiff Matrix" << std::endl;
    GlobalSystem.Print();

    F = RightHandSide;
    GlobalSystem.Solve_LU(F);
    Solution = F;
}

//void PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const{
//
//}

//void PostProcessError(VecDouble error, std::ostream &out, PostProcess &defPostProc) const {
//    
//}
