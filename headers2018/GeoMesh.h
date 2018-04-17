//
//  GeoMesh.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoMesh_h
#define GeoMesh_h

#include "GeoNode.h"
#include "GeoElement.h"
#include <string>

class GeoMesh
{
    
    
    /// vector of nodes
    std::vector<GeoNode> Nodes;
    
    /// vector of element pointers
    std::vector<GeoElement *> Elements;
    
public:
    
    GeoMesh()
    {
        
    }
    
    GeoMesh(const GeoMesh &);
    
    GeoMesh &operator=(const GeoMesh &);
    
    void SetNumNodes(int nnodes);
    
    void SetNumElements(int numelements);
    
    int NumNodes();
    
    int NumElements();
    
    GeoNode &Node(int node);
    
    void SetElement(int elindex, GeoElement *gel);
    
    GeoElement *Element(int elindex);
    
    void BuildConnectivity();
    
    
    
};
#endif /* GeoMesh_h */
