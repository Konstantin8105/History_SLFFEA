/*
    This is the include file "plconst.h" for the finite element progam 
    which uses plate elements.

                Updated 8/5/00

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	pt5  	.5 
#define	pt25 	.25
#define	pt1667  .166666666667
#define	sq3    1.732050808 

#define nsd         2			     /* spatial dimensions per node */
#define ndof        3			     /* degrees of freedom per node */
#define npel        4			     /* nodes per element */
#define neqel       npel*ndof		     /* degrees of freedom per element */
#define num_int     4			     /* number of integration points */
#define neqlsq      neqel*neqel 	     /* neqel squared */
#define lenb        3*neqel*num_int	     /*  */
#define sdim        5                        /* stress dimensions per element */
#define soB         sdim*neqel               /* size of B matrix */
#define MsoB        (nsd + 1)*neqel          /* size of B_mass matrix */
#define soshb       (nsd + 1)*npel*(num_int) /* size of shl and shg matrix
						for 2X2 PT. Gauss for bending */
#define soshs       (nsd +1)*npel*(1)        /* size of shl and shg matrix
						for 1X1 PT. Gauss for shear */
#define sosh_node2  nsd*2*num_int            /* size of shl_node2 */
#define KB            1024.0                      /* number of bytes in kilobyte */
#define MB            1.0486e+06                  /* number of bytes in megabyte */

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


/*
  title               problem title
  numel               number of elements
  numnp               number of nodal points
  nmat                number of materials
  dof                 total number of degrees of freedom
  coord               *(coord + 0) x coordinate of node
                      *(coord + 1) y coordinate of node
  connect             (0-3) connectivity array
  matl                material structure
  Emod                young's modulus
  nu                  poisson's ratio
  thick               thickness of plate
  shear               shear correction factor
  force               *(force + 0) z component of applied load
                      *(force + 1) phi1 component of applied load
                      *(force + 2) phi2 component of applied load
  analysis_flag       1 calculate unknown displacemnts
                      2 calculate reaction forces
  lin_algebra_flag    0 if numel <= 750 elements, use LU Decomposition for
                        displacements
                      1 if numel > 750 elements, use conjugate
                        gradient method for displacements
*/

