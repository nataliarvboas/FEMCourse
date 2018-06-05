//
//  L2Projection.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef L2Projection_h
#define L2Projection_h

#include "MathStatement.h"
#include "DataTypes.h"
#include  "IntPointData.h"
#include <functional>

class L2Projection : public MathStatement
{
    
    // L2 projection matrix
    Matrix projection;
    
    // Force funtion related to L2 projection math statement
    std::function<void(const VecDouble &co, VecDouble &result)> forceFunction;
    
public:
    
    enum PostProcVar {ENone, ESol, EDSol};
    
    // Default constructor of L2Projection
    L2Projection();
    
    // Constructor of L2Projection
    L2Projection(int materialid, Matrix &perm);
    
    // Copy constructor of L2Projection
    L2Projection(const L2Projection &copy);
    
    // Operator of copy
    L2Projection &operator=(const L2Projection &copy);
    
    // Method for creating a copy of the element
    virtual L2Projection *Clone() const;
    
    // Default contructor of L2Projection
    virtual ~L2Projection();
    
    // Return the L2 projection matrix
    Matrix GetProjectionMatrix() const;
    
    // Set the L2 projection matrix
    void SetProjectionMatrix(const Matrix &proj);

    // Return the force function related to L2 projection math statement
    std::function<void(const VecDouble &co, VecDouble &result)> GetForceFunction() const
    {
        return forceFunction;
    }

    // Set the force function related to L2 projection math statement
    void SetForceFunction(const std::function<void(const VecDouble &co, VecDouble &result)> &f)
    {
        forceFunction = f;
    }
    
    // Return the number of errors
    virtual int NEvalErrors() const;
    
    // Return the number of state variables
    virtual int NState() const;
    
    // Return the variable index associated with the name
    virtual int VariableIndex(const std::string &name);
    
    // Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
    virtual int NSolutionVariables(const PostProcVar var);
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, double weight, Matrix &EK, Matrix &EF) const;
    
    // Method to implement error over element's volume
    virtual void ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const;
    
    // Prepare and print post processing data
    virtual std::vector<double> PostProcessSolution(const IntPointData &integrationpointdata, const PostProcVar var) const;
    
    
};
#endif /* L2Projection_h */