/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

//#include "DataTypes.h"
#include "GeoNode.h"

void GeoNode::Print(std::ostream &out) {


    out << "Coordinates:\t";
    int i;
    out << xco[0];
    for (i = 1; i < xco.size(); i++) {
        out << "\t" << xco[i];
    }
    out << std::endl;
}