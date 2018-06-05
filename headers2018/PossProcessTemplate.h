//
//  PostProcessTemplate.h
//  FemCourse
//
//  Created by Philippe Devloo on 15/05/18.
//

#ifndef PostProcessTemplate_h
#define PostProcessTemplate_h

#include "MathStatement.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "IntPointData.h"
#include "PostProcess.h"
#include <list>

class PostProcess;

template<class math>
class PostProcessTemplate: public PostProcess
{

    public:
    
    PostProcessTemplate() : PostProcess(){
    }
    
    PostProcessTemplate(const PostProcessTemplate &copy) : PostProcess(copy) {
    }
    
    ~PostProcessTemplate(){
        
    }
    
    
    PostProcessTemplate &operator=(const PostProcessTemplate &cp){
        return *this;
    }
    
    PostProcessTemplate(Analysis *Ref) : PostProcess(Ref){
    }
    
    
    std::vector<typename math::PostProcVar> scalarvariables;
    std::vector<typename math::PostProcVar> vectorvariables;
    //TVec<typename math::PostProcVar> teste;
    
    virtual void AppendVariable(typename math::PostProcVar obj){
        math Statement;
        int nsol = Statement.NSolutionVariables(obj);
        
        if (nsol==1) {
            scalarvariables.push_back(obj);
        }else{
            vectorvariables.push_back(obj);
        }
        
    }
    
    std::vector<double> PostProcResult(MathStatement &mathStatement, unsigned int varIndex, const IntPointData &data) const {
        math *locptr = dynamic_cast<math *> (&mathStatement);
        
        if(!locptr) DebugStop();
        const int numScalarVariables = NumScalarVariables();
        return locptr->PostProcessSolution(data, varIndex < numScalarVariables? scalarvariables[varIndex] : vectorvariables[varIndex-numScalarVariables]);
    }
    
    inline unsigned int NumScalarVariables() const {
        return scalarvariables.size();
    }
    
    inline unsigned int NumVectorVariables() const {
        return vectorvariables.size();
    }
    
    inline unsigned int NumVariables() const {
        return NumScalarVariables() + NumVectorVariables();
    }
};
#endif /* PostProcessTemplate_h */

template class PostProcessTemplate<Poisson>;
template class PostProcessTemplate<L2Projection>;