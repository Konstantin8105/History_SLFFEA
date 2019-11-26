/*
    This file contains the structures of the shell FEM code.

	Updated 9/22/06

    SLFFEA source file
    Version:  1.4
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shconst.h"

typedef struct {
	double xx,yy,xy,zx,yz,I,II,III;
} SDIM;

typedef struct {
	double x, y, z, phix, phiy;
} XYZPhiF;

typedef struct {
	int x, y, z, phix, phiy;
} XYZPhiI;

typedef struct {
	XYZPhiI *num_fix;
	int   *num_force;
	XYZPhiI *fix;
	int   *force;
} BOUND;

typedef struct {
	double E;
	double nu;
	double rho;
	double thick;
	double shear;
	int extrathick;
} MATL;

typedef struct {
	double *bend;
	double *shear;
	double *bend_z;
	double *shear_z;
} SH;

typedef struct {
	double hat[nsd*npell];
	double bar[nsd*npell];
} XL;

typedef struct {
	double hat[soxshat];
	double bar[soxsbar];
} XLXS;

typedef struct {
	double *l_shear;
	double *l_bend;
	double *f_shear;
} ROTATE;

typedef struct {
	SDIM pt[num_int];
} STRESS;

typedef struct {
	SDIM pt[num_int];
} STRAIN;
