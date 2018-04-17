//
//  TVKGeoMesh.cpp
//  Fem2018
//
//  Created by Philippe Devloo on 17/04/18.
//

#include <stdio.h>

#include "VTKGeoMesh.h"
#include "GeoMesh.h"
#include "tpanic.h"
#include <fstream>
#include <sstream>

static int GetVTK_ElType(ElementType ElType)
{
    
    
    int elType = -1;
    switch (ElType)
    {
        case(EPoint):
        {
            elType = 1;
            break;
        }
        case(EOned):
        {
            elType = 3;
            break;
        }
        case (ETriangle):
        {
            elType = 5;
            break;
        }
        case (EQuadrilateral):
        {
            elType = 9;
            break;
        }
        case (ETetraedro):
        {
            elType = 10;
            break;
        }
        case (EPiramide):
        {
            elType = 14;
            break;
        }
        case (EPrisma):
        {
            elType = 13;
            break;
        }
        case (ECube):
        {
            elType = 12;
            break;
        }
        default:
        {
            std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
            break;
        }
    }
    if(elType == -1)
    {
        std::cout << "Element type not found on " << __PRETTY_FUNCTION__ << std::endl;
        std::cout << "MIGHT BE CURVED ELEMENT (quadratic or quarter point)" << std::endl;
        DebugStop();
    }
    
    return elType;
}

/**
 * Generate an output of all geomesh to VTK
 */
void VTKGeoMesh::PrintGMeshVTK(GeoMesh * gmesh, std::string &filename)
{
    std::ofstream file(filename);
    file.clear();
    int64_t nelements = gmesh->NumElements();
    GeoElement *gel;
    
    std::stringstream node, connectivity, type, material;
    
    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;
    
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";
    
    int64_t actualNode = -1, size = 0, nVALIDelements = 0;
    
    for(int64_t el = 0; el < nelements; el++)
    {
        gel = gmesh->Element(el);
        if(!gel )//|| (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
       
        ElementType elt = gel->Type();
        int elNnodes = gel->NCornerNodes();
        
        size += (1+elNnodes);
        connectivity << elNnodes;
        
        for(int t = 0; t < elNnodes; t++)
        {
            for(int c = 0; c < 3; c++)
            {
                int nodeindex = gel->NodeIndex(t);
                double coord = gmesh->Node(nodeindex).Coord(c);
                node << coord << " ";
            }
            node << std::endl;
            
            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;
        
        int elType = GetVTK_ElType(gel->Type());
        type << elType << std::endl;
        
        material << gel->Material() << std::endl;
        
        nVALIDelements++;
    }
    node << std::endl;
    actualNode++;
    file << actualNode << " float" << std::endl << node.str();
    
    file << "CELLS " << nVALIDelements << " ";
    
    file << size << std::endl;
    file << connectivity.str() << std::endl;
    
    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << type.str() << std::endl;
    
    file << "CELL_DATA" << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 1" << std::endl;
    file << "material 1 " << nVALIDelements << " int" << std::endl;
    file << material.str();
    
    file.close();
}


