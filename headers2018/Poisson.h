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
    
    enum PostProcVar {ENone, ESol, EDSol, EFlux, EForce, ESolExact, EDSolExact};
    
    // Default constructor of Poisson
    Poisson();
    
    // Constructor of Poisson
    Poisson(int materialid, Matrix &perm);
    
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
    
    // returns the integrable dimension of the material
    int Dimension() const {return 2;}
    
    virtual int NEvalErrors() const;
    
    // Return the number of state variables
    virtual int NState() const;
    
    // Return the variable index associated with the name
    virtual int VariableIndex(const std::string &name);
    
    // Return the number of variables associated with the variable indexed by var. Param var Index variable into the solution, is obtained by calling VariableIndex
    virtual int NSolutionVariables(const PostProcVar var);
    
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, double weight , Matrix &EK, Matrix &EF) const;
    
    // Method to implement error over element's volume
    virtual void ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const;
    
    
    // Prepare and print post processing data
    virtual std::vector<double> PostProcessSolution(const IntPointData &integrationpointdata, const PostProcVar var) const;


};
#endif /* Poisson_h */