/*
    This file contains the structures of the plate FEM code.

	Updated 5/25/00

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
#include "plconst.h"

typedef struct {
        double xx,yy,xy,I,II;
} MDIM;

typedef struct {
        double zx,yz;
} SDIM;

typedef struct {
        double z, phix, phiy;
} ZPhiF;

typedef struct {
        int z, phix, phiy;
} ZPhiI;

typedef struct {
	ZPhiI *num_fix;
	int   *num_force;
	ZPhiI *fix;
	int   *force;
} BOUND;

typedef struct {
        double E;
        double nu;
        double rho;
        double thick;
        double shear;
} MATL;

typedef struct {
        double *bend;
	double *shear;
} SH;

typedef struct {
        MDIM pt[num_int];
} MOMENT;

typedef struct {
        SDIM pt[num_int];
} STRESS;

typedef struct {
        MDIM pt[num_int];
} CURVATURE;

typedef struct {
        SDIM pt[num_int];
} STRAIN;
