/*
    This file contains the structures of the graphics program
    for truss elements.

	Updated 9/22/00

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
#include "../truss/tsconst.h"

typedef struct {
        int xx;
} ISDIM;

typedef struct {
        int xx;
} ISTRESS;

typedef struct {
        int xx;
} ISTRAIN;

/* The structure below is a repeat of XYZPhiF found in ../truss/tsstruct.h.
   I cannot simply include tsstruct.h in here because tsstruct.h is
   already included in other modules which tsstrcgr.h is included in
   and this causes a redundancy which is not allowed. */

typedef struct {
        double x, y, z;
} XYZF_GR;

