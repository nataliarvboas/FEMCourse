//
//  PostProcess.h
//  FemCourse
//
//  Created by Philippe Devloo on 15/05/18.
//
#include <functional>
#include "DataTypes.h"

#ifndef PostProcess_h
#define PostProcess_h

class Analysis;

class PostProcess
{
protected:
    
    Analysis *Reference;

    // Pointer to Exact solution function, it is necessary to calculating errors
    std::function<void (const VecDouble &loc, VecDouble &result, Matrix &deriv)> fExact;
    
public:
    
    PostProcess(){
        Reference=0;
    }

    PostProcess(const PostProcess &copy){
        Reference=copy.Reference;
    }
    
    ~PostProcess(){
        
    }
    
    
    PostProcess &operator=(const PostProcess &cp){
        Reference=cp.Reference;
        return *this;
    }
    
    PostProcess(Analysis *Ref){
        Reference=Ref;
    }
    
    virtual void Write(std::string filename){
        
    }
    
    virtual void SetExact(std::function<void (const VecDouble &loc, VecDouble &result, Matrix &deriv)> Exact){
        fExact=Exact;
    }
    
    virtual std::function<void (const VecDouble &loc, VecDouble &result, Matrix &deriv)> GetExact(){
        return fExact;
    }
    
    
};
#endif /* PostProcess_h */