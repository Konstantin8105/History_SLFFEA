/*
     This library function assembles the B matrix in [Btrans][D][B]
     for triangle elements.

     It is a C translation of the subroutine QDCB from the book
     "The Finite Element Method" by Thomas Hughes, page 780.

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

int triangleB(double *shg,double *B)
{
/*
 ....  SET UP THE STRAIN-DISPLACEMENT MATRIX "B" FOR
       TWO-DIMENSIONAL CONTINUUM ELEMENTS
       FOR THE D MATRIX
		Updated 4/20/99
*/
	int i,i2,i2m1;
	for( i = 0; i < npel; ++i )
	{
        	i2      =ndof*i+1;
        	i2m1    =i2-1;

        	*(B+i2m1)         = *(shg+i);
             /* *(B+i2)           = 0.0;
                *(B+neqel*1+i2m1) = 0.0; */
        	*(B+neqel*1+i2)   = *(shg+npel*1+i); 
        	*(B+neqel*2+i2m1) = *(shg+npel*1+i);
        	*(B+neqel*2+i2)   = *(shg+i);
	}
        return 1;
}

int triangleB_mass(double *shg,double *B_mass)
{
/*
     This library function assembles the B_mass matrix in
     [B_mass trans][B_mass] for triangle elements.

                Updated 8/17/00
*/

        int i,i2,i2m1,i2m2;
        for( i = 0; i < npel; ++i )
        {
                i2      =ndof*i+1;
                i2m1    =i2-1;

                *(B_mass+i2m1)       = *(shg+i);
                *(B_mass+neqel*1+i2) = *(shg+i);
        }
        return 1;
}
