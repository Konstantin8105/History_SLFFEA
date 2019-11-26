/* This program reapplies the boundary conditions for the
   displacement when the conjugate gradient method is used.
   This is for a finite element program which does analysis
   on a brick element.

                        Updated 5/4/00 

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
#include "brconst.h"
#if BRICK1
#include "brstruct.h"
#endif
#if BRICK2
#include "../brick2/br2struct.h"
#endif

int brBoundary ( double *U, BOUND bc )
{
	int i,dum;

        for( i = 0; i < bc.num_fix[0].x; ++i )
        {
            dum = bc.fix[i].x;
            *(U+ndof*dum) = 0.0;
        }
        for( i = 0; i < bc.num_fix[0].y; ++i )
        {
            dum = bc.fix[i].y;
            *(U+ndof*dum+1) = 0.0;
        }
        for( i = 0; i < bc.num_fix[0].z; ++i )
        {
            dum = bc.fix[i].z;
            *(U+ndof*dum+2) = 0.0;
        }
	return 1;
}