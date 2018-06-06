/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "CompElementTemplate.h"
#include "CompElement.h"
#include "CompMesh.h"
#include "DataTypes.h"
#include "CompMesh.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTetrahedron.h"
#include "ShapeTriangle.h"
#include "GeoElement.h"
#include "MathStatement.h"
#include "GeoElementSide.h"

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate() : dofindexes(0) {
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh, GeoElement *geo) : CompElement(ind, cmesh, geo) {
    int64_t nel = 0;
    nel = cmesh->GetElementVec().size();
    cmesh->SetNumberElement(nel);
    cmesh->SetElement(ind, this);
    intrule.SetOrder(2 * cmesh->GetDefaultOrder());
    this->SetIntRule(&intrule);
    this->SetIndex(ind);
    this->SetCompMesh(cmesh);
    this->SetGeoElement(geo);
    geo->SetReference(this);

    int matid = geo->Material();
    MathStatement *mat = cmesh->GetMath(matid);
    this->SetStatement(mat);

    int nsides = geo->NSides();

    for (int i = 0; i < nsides; i++) {
        GeoElementSide gelside(this->GetGeoElement(), i);
        GeoElementSide neighbour = gelside.Neighbour();
        this->SetNDOF(nsides);

        if (gelside != neighbour && neighbour.Element()->GetIndex() < ind) {
            CompElement *cel = neighbour.Element()->GetReference();
            int64_t dofindex = cel->GetDOFIndex(neighbour.Side());
            this->SetDOFIndex(i, dofindex);
        } else {
            DOF dof;
            int order = cmesh->GetDefaultOrder();
            int nshape = this->ComputeNShapeFunctions(i, order);
            int nstate = this->GetStatement()->NState();
            dof.SetNShapeStateOrder(nshape, nstate, order);

            int64_t ndof = cmesh->GetNumberDOF();
            ndof++;
            cmesh->SetNumberDOF(ndof);

            int64_t dofindex = ndof - 1;
            this->SetDOFIndex(i, dofindex);

            cmesh->SetDOF(ndof - 1, dof);
        }
    }
}

template<class Shape >
CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate & copy) {
    dofindexes = copy.dofindexes;
    intrule = copy.intrule;
}

template<class Shape >
CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate & copy) {
    dofindexes = copy.dofindexes;
    intrule = copy.intrule;
    return *this;
}

template<class Shape >
CompElementTemplate<Shape>::~CompElementTemplate() {
}

template<class Shape >
CompElement * CompElementTemplate<Shape>::Clone() const {
    return new CompElementTemplate(*this);
}

template<class Shape>
void CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix & dphi) const {
    int order = 1;
    int nsides = this->GetGeoElement()->NCornerNodes();
    VecInt orders(nsides, order);

    Shape::Shape(intpoint, orders, phi, dphi);
}

template<class Shape>
void CompElementTemplate<Shape>::GetMultiplyingCoeficients(VecDouble & coefs) const {
    



    
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions() const {
    int order = 1;
    int nsides = this->GetGeoElement()->NSides();
    VecInt orders(nsides, order);

    return Shape::NShapeFunctions(orders);
}

template<class Shape>
void CompElementTemplate<Shape>::SetNDOF(int64_t ndof) {
    dofindexes.resize(ndof);
}

template<class Shape>
void CompElementTemplate<Shape>::SetDOFIndex(int i, int64_t dofindex) {
    dofindexes[i] = dofindex;
}

template<class Shape >
int64_t CompElementTemplate<Shape>::GetDOFIndex(int i) {
    return dofindexes[i];
}

template<class Shape>
int CompElementTemplate<Shape>::NDOF() const {
    return dofindexes.size();
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) const {
    CompMesh *compmesh = GetCompMesh();
    return compmesh->GetDOF(doflocindex).GetNShape();
}

template<class Shape>
int CompElementTemplate<Shape>::ComputeNShapeFunctions(int doflocindex, int order) {
    dofindexes.resize(doflocindex + 1);
    dofindexes[doflocindex] = doflocindex;
    return Shape::NShapeFunctions(doflocindex, order);

}

template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;