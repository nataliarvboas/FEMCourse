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


template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate() {
}

template<class Shape>
CompElementTemplate<Shape>::CompElementTemplate(int64_t ind, CompMesh *cmesh, GeoElement *geo) : CompElement(ind, cmesh, geo) {
    int nel = cmesh->GetElementVec().size();
    cmesh->GetElementVec().resize(nel+1);
    cmesh->GetElementVec()[nel] = this;
    this->SetIndex(nel);  
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
    //return new CompElementTemplate(*this);
}

template<class Shape>
void CompElementTemplate<Shape>::ShapeFunctions(const VecDouble &intpoint, VecDouble &phi, Matrix &dphi) {
    int order = intrule.GetOrder();    
    VecInt orders(1, order);
    
    Shape::Shape(intpoint, orders, phi, dphi);    
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions() {
    int order = intrule.GetOrder();
    VecInt orders(1, order);
    
    return Shape::NShapeFunctions(orders);
}

template<class Shape>
int CompElementTemplate<Shape>::NDOF() {
    return dofindexes.size();
}

template<class Shape>
int CompElementTemplate<Shape>::NShapeFunctions(int doflocindex) {
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