//
//  Poisson.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef Poisson_h
#define Poisson_h

#include "MathStatement.h"
#include "DataTypes.h"
#include  "IntPointData.h"
#include <functional>

class Poisson : public MathStatement
{
    // Permeability matrix
    Matrix permeability;
    
    // Force funtion related to Poisson math statement
    std::function<void(const VecDouble &co, VecDouble &result)> forceFunction;
    
public:
    
    enum PostProcVar {ENone, ESol, EDSol, EFlux, EForce};
    
    // Default constructor of Poisson
    Poisson();
    
    // Constructor of Poisson
    Poisson(Matrix &perm);
    
    // Copy constructor of Poisson
    Poisson(const Poisson &copy);
    
    // Operator of copy
    Poisson &operator=(const Poisson &copy);
    
    // Method for creating a copy of the element
    virtual Poisson *Clone() const;

    // Destructor of Poisson
    virtual ~Poisson();
    
    // Return the permeability matrix
    Matrix GetPermeability() const;

    // Set the permeability matrix
    void SetPermeability(const Matrix &perm);
    
    // Return the force function related to Poisson math statement
    std::function<void(const VecDouble &co, VecDouble &result)> GetForceFunction() const
    {
        return forceFunction;
    }
    
    // Set the force function related to Poisson math statement
    void SetForceFunction(const std::function<void(const VecDouble &co, VecDouble &result)> &f)
    {
        forceFunction = f;
    }
    
    // Return the number of state variables
    virtual int NState() const;
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, double weight , Matrix &EK, Matrix &EF) const;
    
    // Prepare and print post processing data
    virtual std::vector<double> PostProcess(const IntPointData &integrationpointdata, const PostProcVar var) const;

};
#endif /* Poisson_h */