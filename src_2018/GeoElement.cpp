/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GeoElement.h"
#include "Geom1d.h"
#include "GeomQuad.h"
#include "GeomTetrahedron.h"
#include "GeomTriangle.h"
#include "tpanic.h"
#include "CompElementTemplate.h"
#include "tpanic.h"
#include "Shape1d.h"
#include "ShapeQuad.h"
#include "ShapeTriangle.h"
#include "ShapeTetrahedron.h"

GeoElement::GeoElement() : GMesh(0), MaterialId(-1), Reference(0), Index(-1) {

}

GeoElement::GeoElement(int materialid, GeoMesh *mesh, int index) {
    MaterialId = materialid;
    GMesh = mesh;
    Index = index;
}

GeoElement::GeoElement(const GeoElement &copy) {
    GMesh = copy.GMesh;
    MaterialId = copy.MaterialId;
    Reference = copy.Reference;
    Index = copy.Index;
}

GeoElement::~GeoElement() {

}

CompElement *GeoElement::CreateCompEl(CompMesh *mesh, int64_t index) {
    GeoElement *gel = this;

    switch (gel->Type()) {
        case EOned:
            return new CompElementTemplate<Shape1d> (index, mesh, gel);
            break;
        case EQuadrilateral:
            return new CompElementTemplate<ShapeQuad> (index, mesh, gel);
            break;
        case ETriangle:
            return new CompElementTemplate<ShapeTriangle> (index, mesh, gel);
            break;
        case ETetraedro:
            return new CompElementTemplate<ShapeTetrahedron> (index, mesh, gel);
    }

}

void GeoElement::Print(std::ostream &out) {
    out << "Number of nodes:\t" << this->NNodes() << std::endl;
    out << "Corner nodes:\t\t" << this->NCornerNodes() << std::endl;
    out << "Nodes indexes\t\t";

    for (int i = 0; i < this->NNodes(); i++) out << this->NodeIndex(i) << "   ";
    out << "\nNumber of sides:\t" << this->NSides() << std::endl;
    out << "Material Id:\t\t" << this->Material() << std::endl;

    int nsides = this->NSides();

    for (int i = 0; i < nsides; i++) {
        out << "Neighbours for side " << i << ":\t";
        GeoElementSide neighbour = Neighbour(i);
        GeoElementSide thisside(this, i);
        if (!(neighbour.Element() != 0 && neighbour.Side() > -1)) {
            out << "No neighbour\n";
        } else {
            while (neighbour.Element() != thisside.Element() || neighbour.Side() != thisside.Side()) {
                out << neighbour.Element()->GetIndex() << "/" << neighbour.Side() << ' ';

                neighbour = neighbour.Neighbour();
            }
            out << std::endl;
        }
    }
}