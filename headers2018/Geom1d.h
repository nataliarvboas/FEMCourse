//
//  Geom1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#ifndef Geom1d_h
#define Geom1d_h

#include "Topology1d.h"

class Geom1d : public Topology1d
{
public:
    
    const int NNodes = 2;
    
    /// Constructor
    Geom1d();
    
    /// destructor
    ~Geom1d();
    
    /// copy constructor
    Geom1d(const Geom1d &copy);
    
    /// operator=
    Geom1d &operator=(const Geom1d &copy);
    
    /// Computes the shape functions associated with the geometric map
    static void Shape(const VecDouble &xi, VecDouble &phi, Matrix &dphi);
    
    /// Computes the value of x for a given point in parameter space as a function of corner coordinates
    static void X(const VecDouble &xi, Matrix &NodeCo, VecDouble &x);
    
    /// Computes the value of x and gradx for a given point in parameter space
    static void GradX(const VecDouble &xi, Matrix &NodeCo, VecDouble &x, Matrix &gradx);
    
    /// Set the node indices of the element
    void SetNodes(const VecInt &nodes);
    
    /// Set the node indices of the element
    void GetNodes(VecInt &nodes);
    
    /// Return the index of a node
    int NodeIndex(int node);
    
protected:
    
    /// indexes of the nodes of the geometry
	VecInt fNodeIndices;
};

#endif /* Geom1d_h */
