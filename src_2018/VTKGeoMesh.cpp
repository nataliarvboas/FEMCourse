//
//  TVKGeoMesh.cpp
//  Fem2018
//
//  Created by Philippe Devloo on 17/04/18.
//

#include <stdio.h>

#include "VTKGeoMesh.h"
#include "GeoMesh.h"
#include "GeoElement.h"
#include "CompElement.h"
#include "CompMesh.h"
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
            elType = 21;
            break;
        }
        case (ETriangle):
        {
            // quadratic triangle
            elType = 22;
            break;
        }
        case (EQuadrilateral):
        {
            // quadratic quad
            elType = 23;
            break;
        }
        case (ETetraedro):
        {
            // quadratic tetrahedra
            elType = 24;
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

static TMatrix NodeCoordinates(ElementType eltype)
{
    double quadco[8][2] = {
        {-1,-1},
        {1,-1},
        {1,1},
        {-1,1},
        {0,-1},
        {1,0},
        {0,1},
        {-1,0}
    };
    double triangle[6][2] = {
        {0,0},
        {1,0},
        {0,1},
        {0.5,0},
        {0.5,0.5},
        {0,0.5}
    };
    double tetra[10][3] = {
        {0,0,0},
        {1,0,0},
        {0,1,0},
        {0,0,1},
        {0.5,0,0},
        {0.5,0.5,0},
        {0,0.5,0},
        {0,0,0.5},
        {0.5,0,0.5},
        {0,0.5,0.5}
    };
    TMatrix result;
    switch (eltype) {
        case EOned:
            result.Resize(3, 1);
            result(0,0) = -1.;
            result(1,0) = 1.;
            result(2,0) = 0.;
            break;
        case EQuadrilateral:
            result.Resize(8, 2);
            for (int i=0; i<8; i++) {
                result(i,0) = quadco[i][0];
                result(i,1) = quadco[i][1];
            }
            break;
        case ETriangle:
            result.Resize(6, 2);
            for (int i=0; i<6; i++) {
                result(i,0) = triangle[i][0];
                result(i,1) = triangle[i][1];
            }
            break;
        case ETetraedro:
            result.Resize(10, 3);
            for (int i=0; i<10; i++) {
                result(i,0) = tetra[i][0];
                result(i,1) = tetra[i][1];
                result(i,2) = tetra[i][2];
            }
            break;
        default:
            break;
    }
    return result;
}
/**
 * Generate an output of all geomesh to VTK
 */
void VTKGeoMesh::PrintGMeshVTK(GeoMesh * gmesh, const std::string &filename)
{
    std::ofstream file(filename);
    file.clear();
    int64_t nelements = gmesh->NumElements();
    GeoElement *gel;
    
    std::stringstream node, connectivity, type, material, elindex;
    
    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "GeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;
    
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";
    
    int64_t actualNode = -1, Size = 0, nVALIDelements = 0;
    
    for(int64_t el = 0; el < nelements; el++)
    {
        gel = gmesh->Element(el);
        if(!gel )//|| (gel->Type() == EOned && !gel->IsLinearMapping()))//Exclude Arc3D and Ellipse3D
        {
            continue;
        }
       
        TMatrix ParamCo = NodeCoordinates(gel->Type());
        int elNnodes = ParamCo.Rows();
        
        Size += (1+elNnodes);
        connectivity << elNnodes;
        
        for(int t = 0; t < elNnodes; t++)
        {
            VecDouble xi(ParamCo.Cols(),0.), xco(3,0.);
            for(int i=0; i< xi.size(); i++) xi[i] = ParamCo(t,i);
            gel->X(xi, xco);
            for (auto x:xco) {
                node << x << " ";
            }
            node << std::endl;
            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;
        
        int elType = GetVTK_ElType(gel->Type());
        type << elType << std::endl;
        
        material << gel->Material() << std::endl;
        elindex << el << std::endl;
        nVALIDelements++;
    }
    node << std::endl;
    actualNode++;
    file << actualNode << " float" << std::endl << node.str();
    
    file << "CELLS " << nVALIDelements << " ";
    
    file << Size << std::endl;
    file << connectivity.str() << std::endl;
    
    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << type.str() << std::endl;
    
    file << "CELL_DATA" << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 1" << std::endl;
    file << "material 1 " << nVALIDelements << " int" << std::endl;
    file << material.str();
    file << "FIELD FieldData 1" << std::endl;
    file << "elindex 1 " << nVALIDelements << " int" << std::endl;
    file << elindex.str();

    file.close();
}


/// Generate an output file for the solution and its gradient
void VTKGeoMesh::PrintCMeshVTK(CompMesh *cmesh, int dim, const std::string &filename)
{
    std::ofstream file(filename);
    file.clear();
    
    
    //Header
    file << "# vtk DataFile Version 3.0" << std::endl;
    file << "TPZGeoMesh VTK Visualization" << std::endl;
    file << "ASCII" << std::endl << std::endl;
    
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS ";

    int64_t actualNode = -1, Size = 0, nVALIDelements = 0;

    std::stringstream node, connectivity, Type, material, elindex, solution, gradsol;
    int64_t nelements = cmesh->GetElementVec().size();
    GeoElement *gel;
    for(auto cel:cmesh->GetElementVec())
    {
        gel = cel->GetGeoElement();
        
        TMatrix ParamCo = NodeCoordinates(gel->Type());
        int elNnodes = ParamCo.Rows();
        
        Size += (1+elNnodes);
        connectivity << elNnodes;
        
        for(int t = 0; t < elNnodes; t++)
        {
            VecDouble xi(ParamCo.Cols(),0.), xco(3,0.);
            for(int i=0; i< xi.size(); i++) xi[i] = ParamCo(t,i);
            gel->X(xi, xco);
            for (auto x:xco) {
                node << x << " ";
            }
            node << std::endl;
            VecDouble sol(1);
            TMatrix dsol(2,1);
            cel->Solution(xi, sol, dsol);
            solution << sol[0] << " " << std::endl;
            int i;
            for (i=0; i<dsol.Rows(); i++) {
                gradsol << dsol(i,0) << " ";
            }
            for( ; i<3; i++) gradsol << "0 ";
            gradsol << std::endl;
            actualNode++;
            connectivity << " " << actualNode;
        }
        connectivity << std::endl;
        solution << std::endl;
        
        int elType = GetVTK_ElType(gel->Type());
        Type << elType << std::endl;
        
        material << gel->Material() << std::endl;
        elindex << cel->GetIndex() << std::endl;
        nVALIDelements++;
    }
    node << std::endl;
    actualNode++;
    file << actualNode << " float" << std::endl << node.str();
    
    file << "CELLS " << nVALIDelements << " ";
    
    file << Size << std::endl;
    file << connectivity.str() << std::endl;
    
    file << "CELL_TYPES " << nVALIDelements << std::endl;
    file << Type.str() << std::endl;
    
    file << "CELL_DATA" << " " << nVALIDelements << std::endl;
    file << "FIELD FieldData 1" << std::endl;
    file << "material 1 " << nVALIDelements << " int" << std::endl;
    file << material.str();
    file << "FIELD FieldData 1" << std::endl;
    file << "elindex 1 " << nVALIDelements << " int" << std::endl;
    file << elindex.str();
    
    (file) << "POINT_DATA " << actualNode << std::endl;
    (file) << "SCALARS " << "Solution" << " float" << std::endl << "LOOKUP_TABLE default\n";
    file << solution.str();
    
    (file) << "VECTORS " << "GradSolution" << " float" << std::endl;
    file << gradsol.str();
    file.close();


    
}