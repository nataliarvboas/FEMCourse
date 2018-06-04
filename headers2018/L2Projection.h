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
#include  <functional>

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
    L2Projection(Matrix &perm);
    
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
    
    // Return the number of state variables
    virtual int NState() const;
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, double weight, Matrix &EK, Matrix &EF) const;
    
    // Prepare and print post processing data
    virtual std::vector<double> PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const;

    
};
#endif /* L2Projection_h */