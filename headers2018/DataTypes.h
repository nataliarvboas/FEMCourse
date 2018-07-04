//
//  IntRule.h
//  FemSC
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef __DATATYPES__
#define __DATATYPES__

#include <cmath>
#include <stdio.h>
#include <vector>
#include "TMatrix.h"
//#include "TPZSSpStructMatrix.h"





typedef TMatrix Matrix;
typedef std::vector<int> VecInt;
typedef std::vector<double> VecDouble;


enum ElementType
{

    /*0*/    fEPoint,
    /*1*/    fEOned,
    /*2*/    fETriangle,
    /*3*/    fEQuadrilateral,
    /*4*/    fETetraedro,
    /*5*/    fEPiramide,
    /*6*/    fEPrisma,
    /*7*/    fECube

         
};


#endif /* defined(__FemSC__DATATYPES__) */
