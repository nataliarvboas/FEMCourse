/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "CompElement.h"
#include "GeoElement.h"
#include "MathStatement.h"
#include "CompElementTemplate.h"
#include "IntPointData.h"
#include "MathStatement.h"
#include "tpanic.h"
#include <math.h> 
using namespace std;

CompElement::CompElement() {
}

CompElement::CompElement(int64_t ind, CompMesh *cmesh, GeoElement *geo) {
    compmesh = cmesh;
    index = ind;
    geoel = geo;
}

CompElement::CompElement(const CompElement &copy) {
    compmesh = copy.compmesh;
    index = copy.index;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
}

CompElement &CompElement::operator=(const CompElement &copy) {
    compmesh = copy.compmesh;
    index = copy.index;
    geoel = copy.geoel;
    intrule = copy.intrule;
    mat = copy.mat;
    return *this;
}

CompElement::~CompElement() {
}

CompElement *CompElement::Clone() const {
    //return this->Clone();
}

MathStatement *CompElement::GetStatement() const {
    return mat;
}

void CompElement::SetStatement(MathStatement *statement) {
    mat = statement;
}

IntRule *CompElement::GetIntRule() const {
    return intrule;
}

void CompElement::SetIntRule(IntRule *irule) {
    intrule = irule;
}

void CompElement::SetIndex(int64_t ind) {
    index = ind;

}

GeoElement *CompElement::GetGeoElement() const {
    return geoel;
}

void CompElement::SetGeoElement(GeoElement *element) {
    geoel = element;
}

CompMesh *CompElement::GetCompMesh() const {
    return compmesh;
}

void CompElement::SetCompMesh(CompMesh *mesh) {
    compmesh = mesh;
}

void CompElement::InitializeIntPointData(IntPointData &data) const {

    const int dim = this->Dimension();
    const int nshape = this->NShapeFunctions();
    const int nstate = this->GetStatement()->NState();

    data.weight = 0;
    data.detjac = 0;
    data.phi.resize(nshape, 0);
    data.dphidx.Resize(dim, nshape);
    data.dphidksi.Resize(dim, nshape);
    data.x.resize(3);
    data.ksi.resize(3, 0);
    data.gradx.Resize(dim, nshape);
    data.axes.Resize(dim, 3);
    data.solution.resize(1, 0);
    data.dsoldksi.Resize(1, 1);
    data.dsoldx.Resize(1, 1);

}

void CompElement::ComputeRequiredData(IntPointData &data, VecDouble &intpoint) const {
    GeoElement *ref;
    CompElement *cel = ref->GetReference();    
    	if (!cel){
		std::cout << "Reference() == NULL\n" << std::endl;
		return;
	}

    ref->GradX(data.ksi, data.x, data.gradx);

    data.detjac = 0.0;
    int nrows = data.gradx.Rows();
    int ncols = data.gradx.Cols();
    int dim = this->Dimension();
    Matrix jac(dim, dim, 0);
    Matrix jacinv(dim, dim, 0);

    switch (dim) {
        case 1:
        {
            data.axes.Resize(dim, 3);
            std::vector<double> v_1(3, 0.);

            for (int i = 0; i < nrows; i++) {
                v_1[i] = data.gradx.GetVal(i, 0);
            }

            double norm_v_1 = 0.;
            for (int i = 0; i < nrows; i++) {
                norm_v_1 += v_1[i] * v_1[i];
            }

            norm_v_1 = sqrt(norm_v_1);
            jac(0, 0) = norm_v_1;
            data.detjac = norm_v_1;
            jacinv(0, 0) = 1.0 / data.detjac;

            data.detjac = fabs(data.detjac);

            for (int i = 0; i < 3; i++) {
                data.axes(0, i) = v_1[i] / norm_v_1;
            }
        }
            break;
        case 2:
        {
            data.axes.Resize(dim, 3);
            std::vector<double> v_1(3, 0.), v_2(3, 0.);
            std::vector<double> v_1_til(3, 0.), v_2_til(3, 0.);

            for (int i = 0; i < nrows; i++) {
                v_1[i] = data.gradx.GetVal(i, 0);
                v_2[i] = data.gradx.GetVal(i, 1);
            }

            double norm_v_1_til = 0.0;
            double norm_v_2_til = 0.0;
            double v_1_dot_v_2 = 0.0;

            for (int i = 0; i < 3; i++) {
                norm_v_1_til += v_1[i] * v_1[i];
                v_1_dot_v_2 += v_1[i] * v_2[i];
            }
            norm_v_1_til = sqrt(norm_v_1_til);

            for (int i = 0; i < 3; i++) {
                v_1_til[i] = v_1[i] / norm_v_1_til;
                v_2_til[i] = v_2[i] - v_1_dot_v_2 * v_1_til[i] / norm_v_1_til;
                norm_v_2_til += v_2_til[i] * v_2_til[i];
            }
            norm_v_2_til = sqrt(norm_v_2_til);


            jac(0, 0) = norm_v_1_til;
            jac(0, 1) = v_1_dot_v_2 / norm_v_1_til;
            jac(1, 1) = norm_v_2_til;

            data.detjac = jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1);

            jacinv(0, 0) = +jac(1, 1) / data.detjac;
            jacinv(1, 1) = +jac(0, 0) / data.detjac;
            jacinv(0, 1) = -jac(0, 1) / data.detjac;
            jacinv(1, 0) = -jac(1, 0) / data.detjac;

            data.detjac = fabs(data.detjac);

            for (int i = 0; i < 3; i++) {
                v_2_til[i] /= norm_v_2_til;
                data.axes(0, i) = v_1_til[i];
                data.axes(1, i) = v_2_til[i];
            }
        }
            break;
        case 3:
        {
            data.axes.Resize(dim, 3);

            for (int i = 0; i < nrows; i++) {
                jac(i, 0) = data.gradx.GetVal(i, 0);
                jac(i, 1) = data.gradx.GetVal(i, 1);
                jac(i, 2) = data.gradx.GetVal(i, 2);
            }

            data.detjac -= jac(0, 2) * jac(1, 1) * jac(2, 0); //- a02 a11 a20
            data.detjac += jac(0, 1) * jac(1, 2) * jac(2, 0); //+ a01 a12 a20
            data.detjac += jac(0, 2) * jac(1, 0) * jac(2, 1); //+ a02 a10 a21
            data.detjac -= jac(0, 0) * jac(1, 2) * jac(2, 1); //- a00 a12 a21
            data.detjac -= jac(0, 1) * jac(1, 0) * jac(2, 2); //- a01 a10 a22
            data.detjac += jac(0, 0) * jac(1, 1) * jac(2, 2); //+ a00 a11 a22

            jacinv(0, 0) = (-jac(1, 2) * jac(2, 1) + jac(1, 1) * jac(2, 2)) / data.detjac; //-a12 a21 + a11 a22
            jacinv(0, 1) = (jac(0, 2) * jac(2, 1) - jac(0, 1) * jac(2, 2)) / data.detjac; //a02 a21 - a01 a22
            jacinv(0, 2) = (-jac(0, 2) * jac(1, 1) + jac(0, 1) * jac(1, 2)) / data.detjac; //-a02 a11 + a01 a12
            jacinv(1, 0) = (jac(1, 2) * jac(2, 0) - jac(1, 0) * jac(2, 2)) / data.detjac; //a12 a20 - a10 a22
            jacinv(1, 1) = (-jac(0, 2) * jac(2, 0) + jac(0, 0) * jac(2, 2)) / data.detjac; //-a02 a20 + a00 a22
            jacinv(1, 2) = (jac(0, 2) * jac(1, 0) - jac(0, 0) * jac(1, 2)) / data.detjac; //a02 a10 - a00 a12
            jacinv(2, 0) = (-jac(1, 1) * jac(2, 0) + jac(1, 0) * jac(2, 1)) / data.detjac; //-a11 a20 + a10 a21
            jacinv(2, 1) = (jac(0, 1) * jac(2, 0) - jac(0, 0) * jac(2, 1)) / data.detjac; //a01 a20 - a00 a21
            jacinv(2, 2) = (-jac(0, 1) * jac(1, 0) + jac(0, 0) * jac(1, 1)) / data.detjac; //-a01 a10 + a00 a11

            data.detjac = fabs(data.detjac);

            data.axes.Zero();
            data.axes(0, 0) = 1.0;
            data.axes(1, 1) = 1.0;
            data.axes(2, 2) = 1.0;
        }
            break;
    }
    this->ShapeFunctions(intpoint, data.phi, data.dphidksi);

    int nshape = this->NShapeFunctions();
    data.dphidx.Resize(dim, nshape);
    int ieq;
    switch (dim) {
        case 0:
        {

        }
            break;
        case 1:
        {
            for (ieq = 0; ieq < nshape; ieq++) {
                data.dphidx(0, ieq) *= jacinv.GetVal(0, 0);
            }
        }
            break;
        case 2:
        {
            for (ieq = 0; ieq < nshape; ieq++) {
                data.dphidx(0, ieq) = jacinv.GetVal(0, 0) * data.dphidksi.GetVal(0, ieq) + jacinv.GetVal(1, 0) * data.dphidksi.GetVal(1, ieq);
                data.dphidx(1, ieq) = jacinv.GetVal(0, 1) * data.dphidksi.GetVal(0, ieq) + jacinv.GetVal(1, 1) * data.dphidksi.GetVal(1, ieq);
            }
        }
            break;
        case 3:
        {
            for (ieq = 0; ieq < nshape; ieq++) {
                data.dphidx(0, ieq) = jacinv.GetVal(0, 0) * data.dphidksi.GetVal(0, ieq) + jacinv.GetVal(1, 0) * data.dphidksi.GetVal(1, ieq) + jacinv.GetVal(2, 0) * data.dphidksi.GetVal(2, ieq);
                data.dphidx(1, ieq) = jacinv.GetVal(0, 1) * data.dphidksi.GetVal(0, ieq) + jacinv.GetVal(1, 1) * data.dphidksi.GetVal(1, ieq) + jacinv.GetVal(2, 1) * data.dphidksi.GetVal(2, ieq);
                data.dphidx(2, ieq) = jacinv.GetVal(0, 2) * data.dphidksi.GetVal(0, ieq) + jacinv.GetVal(1, 2) * data.dphidksi.GetVal(1, ieq) + jacinv.GetVal(2, 2) * data.dphidksi.GetVal(2, ieq);
            }
        }
            break;
    }
    data.x.resize(3, 0.0);
    ref->X(intpoint, data.x);
}

void CompElement::CalcStiff(Matrix &ek, Matrix &ef) const {
    MathStatement *material = this->GetStatement();
    if (!material) {
        std::cout << "Error at CompElement::CalcStiff" << std::endl;
        return;
    }

    IntPointData data;
    this->InitializeIntPointData(data);

    int dim = this->Dimension();
    std::vector<double> intpoint(dim, 0.);
    double weight = 0.;

    IntRule *intrule = this->GetIntRule();

    int intrulepoints = intrule->NPoints();
    for (int int_ind = 0; int_ind < intrulepoints; ++int_ind) {
        intrule->Point(int_ind, intpoint, weight);

        this->ComputeRequiredData(data, intpoint);
        weight *= fabs(data.detjac);

        material->Contribute(data, weight, ek, ef);
    }
}