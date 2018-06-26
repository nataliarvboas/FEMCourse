//
//  MathStatement.h
//  FemCourse
//
//  Created by Philippe Devloo on 24/04/18.
//

#ifndef MathStatement_h
#define MathStatement_h

#include "DataTypes.h"
#include "IntPointData.h"
#include "PostProcess.h"

class MathStatement
{
    
    int MathDim;
    
    // Math statement ID
    int matid = 0;
    
public:
    
//    enum PostProcVar;
    
    static double gBigNumber;
    
    // Constructor of MathStatement
    MathStatement();
    
    // Copy constructor of MathStatement
    MathStatement(const MathStatement &copy);
    
    // Operator of copy
    MathStatement &operator=(const MathStatement &copy);
    
    // Destructor of MathStatement
    virtual ~MathStatement();
    
    // Method for creating a copy of the element
    virtual MathStatement *Clone() const = 0;
    
    // Return the number of state variables
    virtual int NState() const =0;
    
    // Return the number of errors
    virtual int NEvalErrors() const = 0;
    
    // Method to implement integral over element's volume
    virtual void Contribute(IntPointData &integrationpointdata, double weight, Matrix &EK, Matrix &EF) const = 0;
    
    // Method to implement error over element's volume
    virtual void ContributeError(IntPointData &integrationpointdata, VecDouble &u_exact, Matrix &du_exact, VecDouble &errors) const = 0;
    
    virtual void SetMatID(int indexmat){
        matid = indexmat;
    }
    
    virtual int GetMatID(){
        return matid;
    }

    virtual void SetDimension(int dim){
        MathDim = dim;
    };
    
    virtual int Dimension() const{
        return MathDim;
    };
    
    // Prepare and print post processing data
    virtual void PostProcessSolution(const IntPointData &integrationpointdata, const int var, VecDouble &sol) const = 0;
    
    
    virtual void Axes2XYZ(const Matrix &dudaxes, Matrix &dudx, const Matrix &axesv, bool colMajor = true) const;
    
    //Method to print MathStatement
    virtual void Print(std::ostream &out);
    
    
};
#endif /* MathStatement_h */