//
//  Analysis.h
//  FemCourse
//
//  Created by Philippe Devloo on 08/05/18.
//

#ifndef Analysis_h
#define Analysis_h

#include "DataTypes.h"
class CompMesh;
class PostProcess;
#include <string>


class Analysis
{
protected:
    
    CompMesh *cmesh;
    
    Matrix Solution;
    
    Matrix GlobalSystem;
    
    Matrix RightHandSide;
    
public:
    
    Analysis();
    
    Analysis(const Analysis &cp);
    
    Analysis &operator=(const Analysis &cp);
    
    Analysis(CompMesh *cmesh);
    
    void SetMesh(CompMesh *cmesh);
    
    CompMesh *Mesh() const;
    
    void RunSimulation();
    
    void PostProcess(std::string &filename, PostProcess &defPostProc) const;
    
};

#endif /* Analysis_h */