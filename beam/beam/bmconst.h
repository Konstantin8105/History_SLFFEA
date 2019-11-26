/*
    This is the include file "bmconst.h" for the finite element progam 
    which uses 3D beam elements.

                Updated 9/13/00

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	pt5           .5 
#define	pt25          .25
#define	pt1667        .166666666667
#define	sq3          1.732050808 
#define	sq3pt33      1.825741858 
#define	sqpt3        0.547722557 
#define	ONE          0.99999 

#define nsd            3                  /* number of spatial dimensions per node */
#define ndof           6                  /* degrees of freedom per node */
#define npel           2                  /* nodes per element */
#define neqel          npel*ndof          /* degrees of freedom per element */
#define num_int        2                  /* number of integration points */
#define num_int_mass   4                  /* number of integration points for mass */
#define nsdsq          nsd*nsd 	          /* nsd squared */
#define neqlsq         neqel*neqel        /* neqel squared */
#define lenb           3*neqel*num_int    /*  */
#define sdim           4                  /* stress dimensions per element */
#define soB            sdim*neqel         /* size of B matrix */
#define KB             1024.0             /* number of bytes in kilobyte */
#define MB             1.0486e+06         /* number of bytes in megabyte */

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

/*
  title               problem title
  numel               number of elements
  numnp               number of nodal points
  nmat                number of materials
  nmode               number of modes  nmode = 0 means standard analysis
  dof                 total number of degrees of freedom
  coord               *(coord + 0) x coordinate of node
                      *(coord + 1) y coordinate of node
                      *(coord + 2) z coordinate of node
  connect (0-7)       connectivity array
  matl                material structure
  Emod                young's modulus
  nu                  poisson's ratio
  area                area
  Iy                  Rotational Inertia about y axis
  Iz                  Rotational Inertia about z axis
  Ip                  Polar moment of Inertia
  el_type             (1) truss
                      (2) beam
  force               *(force + 0) x component of applied load
                      *(force + 1) y component of applied load
                      *(force + 2) z component of applied load
                      *(force + 3) applied moment about x
                      *(force + 4) applied moment about y
                      *(force + 5) applied moment about z
  analysis_flag       1 calculate unknown displacemnts
                      2 calculate reaction forces
  lin_algebra_flag    0 if numel <= 500 elements, use LU Decomposition for
                        displacements
                      1 if numel > 500 elements, use conjugate
                        gradient method for displacements
*/
