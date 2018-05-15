/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"

L2Projection::L2Projection() {
}

L2Projection::L2Projection(Matrix &perm) {
    projection = perm;
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    return *this;
}

L2Projection *L2Projection::Clone() const {
    //return new L2Projection(*this);
}

L2Projection::~L2Projection() {
}

Matrix L2Projection::GetProjectionMatrix() const {
    return projection;
}

void L2Projection::SetProjectionMatrix(const Matrix &proj) {
    projection = proj;
}

int L2Projection::NState() const {
 
}

void L2Projection::Contribute(IntPointData &integrationpointdata, double weight, Matrix &EK, Matrix &EF) const {
    
}

void L2Projection::PostProcess(IntPointData &integrationpointdata, const std::string &variable, VecDouble &postprocvalue) const {
}