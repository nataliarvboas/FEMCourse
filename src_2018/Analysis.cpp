/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Analysis.h"
#include "Assemble.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "PostProcess.h"
#include "tpanic.h"
#include "MathStatement.h"
#include "VTKGeoMesh.h"
#include "PostProcessTemplate.h"

using namespace std;

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
    Matrix F;
    Assemble assemb(cmesh);
    assemb.Compute(GlobalSystem, RightHandSide);

    //    std::cout << "\nGlobal Stiff Matrix" << std::endl;
    //    GlobalSystem.Print();
    //
    //    std::cout << "\nLoad Vector:" << std::endl;
    //    RightHandSide.Print();

    int nrows = RightHandSide.Rows();
    int ncols = RightHandSide.Cols();
    F.Resize(nrows, ncols);
    F = RightHandSide;

    GlobalSystem.Solve_LU(F);
    Solution = F;

    int solsize = Solution.Rows();
    VecDouble sol(solsize);

    for (int i = 0; i < solsize; i++) {
        sol[i] = Solution(i, 0);
    }
    cmesh->LoadSolution(sol);
}

void Analysis::PostProcessSolution(const std::string &filename, PostProcess &defPostProc) const {
    VTKGeoMesh::PrintSolVTK(cmesh, defPostProc, filename);
}

VecDouble Analysis::PostProcessError(std::ostream &out, PostProcess &defPostProc) const {

    VecDouble values(10, 0.);
    VecDouble errors(10, 0.);
    std::function<void (const VecDouble &loc, VecDouble &result, Matrix & deriv) > fExact;

    int64_t nel = cmesh->GetElementVec().size();

    for (int64_t i = 0; i < nel; i++) {
        CompElement *el = cmesh->GetElement(i);
        if (el) {
            if (el->GetStatement()->GetMatID() == 0) {
                std::fill(errors.begin(), errors.end(), 0.);
                fExact = defPostProc.GetExact();
                el->EvaluateError(fExact, errors);
                int nerrors = errors.size();
                values.resize(nerrors, 0.);
                for (int ier = 0; ier < nerrors; ier++) {
                    values[ier] += errors[ier] * errors[ier];
                }
            }
        }
    }

    int nerrors = errors.size();
    VecDouble ervec(nerrors, -10.0);

    if (nerrors < 3) {
        out << endl << "Analysis::PostProcessError - At least 3 norms are expected." << endl;
        out << endl << "############" << endl;
        for (int ier = 0; ier < nerrors; ier++)
            out << endl << "error " << ier << "  = " << sqrt(values[ier]);
    } else {
        out << "############" << endl;
        out << "Norma L2 de u: " << sqrt(values[0]) << endl;
        out << "Norma L2 de gradu: " << sqrt(values[1]) << endl;
        out << "Norma H1 de u: " << sqrt(values[2]) << endl;
        for (int ier = 3; ier < nerrors; ier++)
            out << "other norms = " << sqrt(values[ier]) << endl;
    }
    // Returns the calculated errors.
    for (int i = 0; i < nerrors; i++) {
        ervec[i] = sqrt(values[i]);
    }
    return ervec;
}

