//
//  GeomTetrahedron.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef GeomTetrahedron_h
#define GeomTetrahedron_h

#include "DataTypes.h"
#include "TopologyTetrahedron.h"

class GeomTetrahedron : public TopologyTetrahedron
{
public:
    
    const static int NNodes = 4;
    
    const static int NSides = 15;
    
    /// Constructor
    GeomTetrahedron();
    
    /// destructor
    ~GeomTetrahedron();
    
    /// copy constructor
    GeomTetrahedron(const GeomTetrahedron &copy);
    
    /// operator=
    GeomTetrahedron &operator=(const GeomTetrahedron &copy);
    
    /// Computes the shape functions associated with the geometric map
    void Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi);
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    void X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x);
    
    /// Computes the value of x and gradx for a given point in parameter space
    void GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx);
    
    /// Set the node indices of the element
    void SetNodes(const VecInt &nodes);
    
    /// Set the node indices of the element
    void GetNodes(VecInt &nodes);
    
    /// Return the index of a node
    int NodeIndex(int node);
    
    /// Return the neighbour along side
    GeoElementSide Neighbour(int side);
    
    /// Initialize the neighbour data structure
    void SetNeighbour(int side, GeoElementSide &neighbour);
    
protected:
    
    /// indexes of the nodes of the geometry
    VecInt fNodeIndices;
    
    /// vector of neighbours
    GeoElementSide fNeighbours[NSides];
};

#endif /* GeomTetrahedron_h */
