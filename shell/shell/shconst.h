/*
    This is the include file "shconst.h" for the finite element progam 
    which uses shell elements.

                Updated 8/9/06

    SLFFEA source file
    Version:  1.4
    Copyright (C) 1999, 2000, 2001, 2002  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/


#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pt5     .5 
#define pt25    .25
#define pt1667  .166666666667
#define sqpt5   .707106781
#define sq3     1.732050808 

#define nsd         3                                    /* spatial dimensions per node */
#define nsdl        2                                    /* isoparametric dimensions in lamina */
#define nsdf        1                                    /* isoparametric dimensions in fiber */
#define ndof        5                                    /* degrees of freedom per node */
#define npelf       2                                    /* nodes per element fiber */
#define npell       4                                    /* nodes per element lamina */
#define npel        npelf*npell                          /* nodes per element */
#define neqel       npell*ndof                           /* degrees of freedom per element */
#define num_ints    2                                    /* number of integ. points on fiber*/
#define num_intb    4                                    /* number of integ. points on lamina */
#define num_int     num_intb*num_ints                    /* number of integ. points on lamina */
#define neqlsq      neqel*neqel                          /* neqel squared */
#define nsdsq       nsd*nsd                              /* nsd squared */
#define lenb        3*neqel*num_intb                     /*  */
#define nremat      5                                    /* number of material properties */
#define nrconnect1  8                                    /* number of possible disp. bound. cond. */
#define nrconnect2  4                                    /* number of possible rotat. bound. cond. */
#define sdim        5                                    /* stress dimensions per element */
#define soB         sdim*neqel                           /* size of B matrix for shell */
#define MsoB        nsd*neqel                            /* size of B_mass matrix */
#define soxsbar     nsd*(nsd-1)                          /* size of xs.bar */
#define soxshat     nsd*nsd                              /* size of xs.hat */

#define soshlb      (nsdl+1)*npell*(num_intb)            /* size of shl matrix
                                                            for 2X2 PT. Gauss for bending */
#define soshls      (nsdl+1)*npell*(1)                   /* size of shl matrix
                                                            for 1X1 PT. Gauss for shear */
#define soshl_zb    (nsdf+1)*npelf*num_ints              /* size of shl_z matrix
                                                            for 2X1 PT. Gauss for bending */
#define soshl_zs    (nsdf+1)*npelf                       /* size of shl_z matrix
                                                            for 1X1 PT. Gauss for shear */
#define soshgb      (nsd+1)*npell*(num_intb)*num_ints    /* size of shg matrix (+1 SRI)
                                                            for 2X2 PT. Gauss for bending */
#define soshgs      (nsd+1)*npell*(1)*num_ints           /* size of shg matrix (+1 SRI)
                                                            for 1X1 PT. Gauss for shear */
#define soshg_zb    (nsd)*npell*(num_intb)*num_ints      /* size of shg_z matrix
                                                            for 2X2 PT. Gauss for bending */
#define soshg_zs    (nsd)*npell*(1)*num_ints             /* size of shg_z matrix
                                                            for 1X1 PT. Gauss for shear */
#define sorlb       nsdsq*num_int                        /* size of lamina rotat. matrix bending */
#define sorls       nsdsq*num_ints                       /* size of lamina rotat. matrix shear */
#define sorfs       nsdsq*npell                          /* size of fiber rotat. matrix shear */
#define zeta        0.00                                 /* reference point on lamina */
#define sosh_node2  nsd*2*num_ints                       /* size of shl_node2 */
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
                      *(coord + 2) z coordinate of node
  connect             (0-3) connectivity array 
  matl                material structure
  Emod                young's modulus
  nu                  poisson's ratio
  shear               shear correction factor
  force               *(force + 0) x component of applied load
                      *(force + 1) y component of applied load
                      *(force + 2) z component of applied load
                      *(force + 3) phi1 component of applied load
                      *(force + 4) phi2 component of applied load
  integ_flag          0 Reduced integration of membrane shear
                      1 Reduced integration of transverse shear
                      2 Reduced integration of membrane and 
                        transverse shear
  analysis_flag       1 calculate unknown displacemnts
                      2 calculate reaction forces
  LU_decomp_flag      0 if numel <= 450 elements, use LU Decomposition for
                        displacements
                      1 if numel > 450 elements, use conjugate
                        gradient method for displacements
*/
