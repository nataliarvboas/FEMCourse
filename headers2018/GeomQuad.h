//
//  GeomQuad.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef GeomQuad_h
#define GeomQuad_h

#include "TopologyQuad.h"
#include "GeoElementSide.h"

class GeomQuad : public TopologyQuad
{
public:
    
    // Constructor
    GeomQuad();
    
    // Destructor
    ~GeomQuad();
    
    // Copy constructor
    GeomQuad(const GeomQuad &copy);
    
    // Operator of copy
    GeomQuad &operator=(const GeomQuad &copy);
    
    // Computes the shape functions associated with the geometric map
    static void Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi);
    
    // Computes the value of x for a given point in parameter space as a function of corner coordinates
    static void X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x);
    
    // Computes the value of x and gradx for a given point in parameter space
    static void GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx);
    
    // Return the number of nodes of the template
    int NumNodes();
    
    // Set the node indices of the element
    void SetNodes(const VecInt &nodes);
    
    // Set the node indices of the element
    void GetNodes(VecInt &nodes);
    
    // Return the index of a node
    int NodeIndex(int node);
    
    // Return the neighbour along side
    GeoElementSide Neighbour(int side);
    
    // Initialize the neighbour data structure
    void SetNeighbour(int side, const GeoElementSide &neighbour);
    
protected:
    
    // Indexes of the nodes of the geometry
    VecInt fNodeIndices;
    
    // Vector of neighbours
    GeoElementSide fNeighbours[nSides];

};

#endif /* GeomQuad_h */