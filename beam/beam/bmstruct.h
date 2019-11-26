/*
    This file contains the structures of the 3-D linear beam
    FEM code.
		Udated 5/22/00

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
#include "bmconst.h"

typedef struct {
        double xx,yy,zz;
} MDIM;

typedef struct {
        double xx;
} SDIM;

typedef struct {
        double x, y, z, phix, phiy, phiz;
} XYZPhiF;

typedef struct {
        double qy, qz;
} QYQZ;

typedef struct {
        int x, y, z, phix, phiy, phiz;
} XYZPhiI;

typedef struct {
	XYZPhiI *num_fix;
        int *num_force;
        int *num_dist_load;
	XYZPhiI *fix;
	int *dist_load;
        int *force;
} BOUND;

typedef struct {
	double E;
	double nu;
	double rho;
	double area;
	double Iy;
	double Iz;
	double Ip;
} MATL;

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

