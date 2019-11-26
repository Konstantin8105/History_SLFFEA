/*
    This file contains the structures of the graphics program
    for triangle elements.

	Updated 8/15/01

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
#include "../tri/trconst.h"

typedef struct {
        int xx,yy,xy,I,II;
} ISTRESS;

typedef struct {
        int xx,yy,xy,I,II;
} ISTRAIN;

/* The structure below is a repeat of XYF found in ../tri/trstruct.h.
   I cannot simply include qdstruct.h in here because qdstruct.h is
   already included in other modules which trstrcgr.h is included in
   and this causes a redundancy which is not allowed. */

typedef struct {
        double x, y;
} XYF_GR;

