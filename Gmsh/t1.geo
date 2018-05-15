/*********************************************************************
 *
 *  Gmsh tutorial 1
 *
 *  Variables, elementary entities (points, lines, surfaces), physical
 *  entities (points, lines, surfaces)
 *
 *********************************************************************/

// The simplest construction in Gmsh's scripting language is the
// `affectation'. The following command defines a new variable `lc':

IsquadQ = 1;

lc = 1e-2;
pr = 1.;
nx= 3;
ny= 3;

Point(1) = {0, 0, 0, lc};

Point(2) = {.1, 0,  0, lc} ;
Point(3) = {.1, .3, 0, lc} ;
Point(4) = {0,  .3, 0, lc} ;

Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Transfinite Line{1,3} = nx Using Progression pr;
Transfinite Line{2,4} = ny Using Progression pr;

Line Loop(1) = {1,2,3,4} ;

Plane Surface(1) = {1} ;

//Verificar se eh quadrilatero ou triangulo aqui:
Recombine Surface {1};


Transfinite Surface{1}={1,2,3,4};

Physical Line("Bottom") = {1} ;
Physical Line("Right") = {2} ;
Physical Line("Top") = {3} ;
Physical Line("Left") = {4} ;
Physical Surface("My surface") = {1} ;

Coherence Mesh;
