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
//#include "TVecNum.h"
#include "TMatrix.h"


typedef TMatrix Matrix;
typedef std::vector<int> VecInt;
typedef std::vector<double> VecDouble;

enum ElementType
{
    /*0*/    EPoint,
    /*1*/    EOned,
    /*2*/    ETriangle,
    /*3*/    EQuadrilateral,
    /*4*/    ETetraedro,
    /*5*/    EPiramide,
    /*6*/    EPrisma,
    /*7*/    ECube
};

#endif /* defined(__FemSC__DATATYPES__) */
