//
//  CompMesh.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef CompMesh_h
#define CompMesh_h

class CompElement;

class MathStatement;

#include "DOF.h"
#include "DataTypes.h"

class CompMesh
{
    // Vector with computational elements objects
    std::vector<CompElement *> compelements;
    
    // Vector with degrees of freedom
    std::vector<DOF> dofs;
    
    // Vector with math statement objects
    std::vector<MathStatement *> mathstatements;
    
public:
    
    // Default constructor of CompMesh
    CompMesh();
    
    // Copy constructor of CompMesh
    CompMesh(const CompMesh &copy);
    
    // Destructor of CompMesh
    ~CompMesh();
    
    // Set the number of computational elements on the grid
    void SetNumberElement(int64_t nelem);
    
    // Set the number of degrees of freedom
    void SetNumberDOF(int64_t ndof);
    
    // Set the number of math statements
    void SetNumberMath(int nmath);
    
    // Set the computational element associated to an index
    void SetElement(int64_t elindex, CompElement *cel);
    
    // Set the degree of freedom associated to an index
    void SetDOF(int64_t index, const DOF &dof);
    
    // Set the math statement object associated to an index
    void SetMathStatement(int index, MathStatement *math);
    
    // Return the degree of freedom index
    DOF &GetDOF(int64_t dofindex);
    
    // Return the computational element associated to an index
    CompElement *GetElement(int64_t elindex) const;
    
    // Return the math statement object associated to an index
    MathStatement *GetMath(int matindex) const;

    // Return the vector with computational elements
    std::vector<CompElement *> GetElementVec() const;
    
    // Return the vector with degrees of freedom
    std::vector<DOF> GetDOFVec() const;
    
    // Return the vector with math statement objects
    std::vector<MathStatement *> GetMathVec() const;
    
    // Set the vector with computational elements
    void SetElementVec(const std::vector<CompElement *> &vec);
    
    // Set the vector with degrees of freedom
    void SetDOFVec(const std::vector<DOF> &dofvec);
    
    // Set the vector with math statement objects
    void SetMathVec(const std::vector<MathStatement *> &mathvec);
    
    // Initialize the datastructure FirstEquation of the DOF objects
    void Resequence();
    
    // Initialize the datastructure FirstEquation of the DOF objects in the order specified by the vector
    void Resequence(VecInt &DOFindices);
    
};

#endif /* CompMesh_h */