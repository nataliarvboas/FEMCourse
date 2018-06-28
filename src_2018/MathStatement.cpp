/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "MathStatement.h"

double MathStatement::gBigNumber = 1.e12;

MathStatement::MathStatement() {
}

MathStatement::MathStatement(const MathStatement &copy) {
    matid = copy.matid;
}

MathStatement &MathStatement::operator=(const MathStatement &copy) {
    matid = copy.matid;
    return *this;
}

MathStatement::~MathStatement() {
}

void MathStatement::Axes2XYZ(const Matrix &dudaxes, Matrix &dudx, const Matrix &axesv, bool colMajor) const {
    Matrix axes(axesv.Rows(), axesv.Cols());
    for (int r = 0; r < axes.Rows(); r++) {
        for (int c = 0; c < axes.Cols(); c++) {
            axes(r, c) = axesv.GetVal(r, c);
        }
    }

    if (colMajor) {
        Matrix axesT;
        axes.Transpose(axesT);

        if (dudx.Rows() != axesT.Rows() || dudx.Cols() != dudaxes.Cols()) {
            dudx.Resize(axesT.Rows(), dudaxes.Cols());
        }
        dudx.Zero();
        axesT.Multiply(dudaxes, dudx, 0);
    } else {
        dudx.Resize(dudaxes.Rows(), axes.Cols());
        dudx.Zero();
        dudaxes.Multiply(axes, dudx, 0);
    }
}

void MathStatement::Print(std::ostream &out) {
    out << "Material: " << this->GetMatID() << std::endl;
    out << "Big number: " << gBigNumber << std::endl << std::endl;

}