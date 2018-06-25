//
//  GeoElement.h
//  FemSC
//
//  Created by Philippe Devloo on 16/04/18.
//

#ifndef GeoElement_h
#define GeoElement_h

#include "DataTypes.h"
#include "CompElement.h"

class CompMesh;

class GeoElementSide;

class GeoMesh;

class GeoElement
{
    
protected:
    // Geometric mesh to which the element belongs
    GeoMesh *GMesh;
    
    // Material ID associated with the element
    int MaterialId;
    
    // Pointer to computational element
    CompElement *Reference;
    
    // Index of the element in the element vector
    int Index;
    
public:

    // Default Constructor of GeoElement
    GeoElement();
    
    // Constructor of GeoElement
    GeoElement(int materialid, GeoMesh *mesh, int index);

    // Copy constructor of GeoElement
    GeoElement(const GeoElement &copy);

    // Destructor of GeoElement
    virtual ~GeoElement();
    
    virtual GeoElement *Clone(GeoMesh *gmesh) const =0;
    
    // Return number of corner nodes
    virtual int NCornerNodes() = 0;
    
    // Return number of nodes
    virtual int NNodes() = 0;
    
    // Return number of sides
    virtual int NSides() = 0;

    // Return the index of an element nodes
    virtual int NodeIndex(int node) = 0;
    
    // Return number fo sides associated with a side
    virtual int NSideNodes(int side) = 0;
    
    // Local node index of a node associated with a side
    virtual int SideNodeIndex(int side, int node) = 0;
    
    /// Set the node indices of the element
    virtual void GetNodes(VecInt &nodes) = 0;
    
    // Return neighbour element of a given side
    virtual GeoElementSide Neighbour(int side) = 0;
    
    // Set neighbour element of a given side
    virtual void SetNeighbour(int side, const GeoElementSide &neigh) = 0;
    
    // Return the enumerated element type
    virtual ElementType Type() = 0;

    //Set reference
    virtual void SetReference(CompElement * elp){
        Reference = elp;
    }
    
    //Get reference
    virtual CompElement *GetReference() const{
        return Reference;
    }
    
    // Set geometric mesh
    void SetMesh(GeoMesh *gmesh)
    {
        GMesh = gmesh;
    }
    
    // Get geometric mesh
    GeoMesh* GetMesh(){
        return GMesh;
    }
    
    // Return material ID
    int Material()
    {
        return MaterialId;
    }
    
    // Create a computational element
    CompElement *CreateCompEl(CompMesh *mesh, int64_t index);
    
    // Set the element index
    void SetIndex(int index)
    {
        Index = index;
    }
    
    // Return the element index
    int GetIndex()
    {
        return Index;
    }
    
    virtual int SideIsUndefined(int side) = 0;
    
    virtual int WhichSide(VecInt &SideNodeIds) = 0;
    
    // Compute x mapping from local parametric coordinates
    virtual void X(const VecDouble &xi,  VecDouble &x) = 0;
    
    // Compute gradient of x mapping from local parametric coordinates
    virtual void GradX(const VecDouble &xi, VecDouble &x, Matrix &gradx) = 0;
    
    // Compute a QR facotrization of the gradient of the mapping function, Q = Jacobian and R = axes
    virtual void Jacobian(const Matrix &gradx, Matrix &jac,Matrix &axes, double &detjac, Matrix &jacinv) = 0;
    
    // Function to print results
    virtual void Print(std::ostream &out);
    

};
#endif /* GeoElement_h */