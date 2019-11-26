/* This program reapplies the boundary conditions for the
   displacement when the conjugate gradient method is used.
   This is for a finite element program which does analysis
   on a wedge element.

                        Updated 3/21/01 

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
#include "weconst.h"
#if WEDGE1
#include "westruct.h"
#endif
#if WEDGE2
#include "../wedge2/we2struct.h"
#endif

int weBoundary ( double *U, BOUND bc )
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
