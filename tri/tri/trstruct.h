/*
    This file contains the structures of the triangle FEM code.

	Updated 8/28/01

    SLFFEA source file
    Version:  1.3
    Copyright (C) 1999, 2000, 2001, 2002  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "trconst.h"

typedef struct {
        double x, y;
} XYF;

typedef struct {
        int x, y;
} XYI;

typedef struct {
	XYI *num_fix;
	int *num_force;
	XYI *fix;
	int *force;
} BOUND;

typedef struct {
        double E;
        double nu;
        double rho;
} MATL;

typedef struct {
        double xx,yy,xy,I,II;
} STRESS;

typedef struct {
        double xx,yy,xy,I,II;
} STRAIN;
