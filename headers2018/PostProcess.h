//
//  PostProcess.h
//  FemCourse
//
//  Created by Philippe Devloo on 15/05/18.
//

#ifndef PostProcess_h
#define PostProcess_h

#include <string>

class Analysis;

class PostProcess
{
    
    Analysis *Reference;
    
    virtual void Write(std::string filename);
        
    
};
#endif /* PostProcess_h */