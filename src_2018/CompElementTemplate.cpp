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

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate() {
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh, GeoElement *geo) : CompElement(ind, cmesh, geo) {
    cmesh->SetGeoMesh(geo->GetMesh());
    int64_t nel = 0;
    nel = cmesh->GetGeoMesh()->NumElements();
    cmesh->SetNumberElement(nel);
    cmesh->SetElement(ind, this);
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(const CompElementTemplate &copy) {
    dofindexes = copy.dofindexes;
    intrule = copy.intrule;
}

template<class Shape>
CompElementTemplate<Shape> &CompElementTemplate<Shape>::operator=(const CompElementTemplate &copy) {
    dofindexes = copy.dofindexes;
    intrule = copy.intrule;
    return *this;
}

template<class Shape>
CompElementTemplate<Shape>::~CompElementTemplate() {
}

template<class Shape>
CompElement *CompElementTemplate<Shape>::Clone() const {
    return new CompElementTemplate(*this);
}

template<class Shape>
void CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) const {
    int order = 1;
    int nsides = this->GetGeoElement()->NCornerNodes();
    VecInt orders(nsides, order);

    Shape::Shape(intpoint, orders, phi, dphi);
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

}

template class CompElementTemplate<Shape1d>;
template class CompElementTemplate<ShapeQuad>;
template class CompElementTemplate<ShapeTriangle>;
template class CompElementTemplate<ShapeTetrahedron>;