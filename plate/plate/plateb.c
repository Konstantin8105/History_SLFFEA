/*
     This library function assembles the B matrix in [Btrans][D][B]
     for plate elements.

     It is based on the subroutine QDCB from the book
     "The Finite Element Method" by Thomas Hughes, page 780.

     SLFFEA source file
     Version:  1.1
     Copyright (C) 1999, 2000  San Le

     The source code contained in this file is released under the
     terms of the GNU Library General Public License.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "plconst.h"

int plateB4pt(double *shg,double *B)
{
/*
 ....  SET UP THE STRAIN-DISPLACEMENT MATRIX "B" FOR
       TWO-DIMENSIONAL CONTINUUM ELEMENTS
       FOR THE D MATRIX FOR 2X2 INTEGRATION

		Updated 5/11/99 
*/
	int i,i2,i2m1,i2m2;
	for( i = 0; i < npel; ++i )
	{
        	i2      =ndof*i+2;
        	i2m1    =i2-1;
        	i2m2    =i2-2;

             /* *(B+i2m2)         = 0.0; 
                *(B+i2m1)         = 0.0; */
        	*(B+i2)           = -*(shg+i);
             /* *(B+neqel*1+i2m2) = 0.0; */
        	*(B+neqel*1+i2m1) = *(shg+npel*1+i); 
             /* *(B+neqel*1+i2)   = 0.0; */
             /* *(B+neqel*2+i2m2) = 0.0; */
        	*(B+neqel*2+i2m1) = *(shg+i);
        	*(B+neqel*2+i2)   = -*(shg+npel*1+i);
	}
        return 1;
}

int plateB1pt(double *shg,double *B)
{
/*
 ....  SET UP THE STRAIN-DISPLACEMENT MATRIX "B" FOR
       TWO-DIMENSIONAL CONTINUUM ELEMENTS
       FOR THE D MATRIX FOR 1 POINT INTEGRATION

		Updated 5/11/99 
*/
	int i,i2,i2m1,i2m2;
	for( i = 0; i < npel; ++i )
	{
        	i2      =ndof*i+2;
        	i2m1    =i2-1;
        	i2m2    =i2-2;

        	*(B+i2m2)         = *(shg+i);
             /* *(B+i2m1)         = 0.0; */
        	*(B+i2)           = *(shg+npel*2+i);
        	*(B+neqel*1+i2m2) = *(shg+npel*1+i);
        	*(B+neqel*1+i2m1) = -*(shg+npel*2+i); 
             /* *(B+neqel*1+i2)   = 0.0; */
	}
        return 1;
}

int plateB4pt_node(double *shg,double *B)
{
/*
 ....  SET UP THE STRAIN-DISPLACEMENT MATRIX "B" FOR
       TWO-DIMENSIONAL CONTINUUM ELEMENTS
       FOR THE D MATRIX FOR 2X2 INTEGRATION.

       THIS DIFFERS FROM "plateB4pt" IN THAT THE CALCULATION
       IS DONE FOR THE ENTIRE "B" MATRIX RATHER THAN JUST
       THE BENDING.  THIS IS USED FOR WHEN THE STRESSES AND
       STRAINS ARE CALCULATED AT THE NODES RATHER THAN GAUSS
       POINTS.

		Updated 6/17/00 
*/
	int i,i2,i2m1,i2m2;
	for( i = 0; i < npel; ++i )
	{
        	i2      =ndof*i+2;
        	i2m1    =i2-1;
        	i2m2    =i2-2;

        	*(B+i2)           = -*(shg+i);
        	*(B+neqel*1+i2m1) = *(shg+npel*1+i); 
        	*(B+neqel*2+i2m1) = *(shg+i);
        	*(B+neqel*2+i2)   = -*(shg+npel*1+i);
        	*(B+neqel*3+i2m2) = *(shg+i);
        	*(B+neqel*3+i2)   = *(shg+npel*2+i);
        	*(B+neqel*4+i2m2) = *(shg+npel*1+i);
        	*(B+neqel*4+i2m1) = -*(shg+npel*2+i); 
	}
        return 1;
}


int plateB_mass(double *shg,double *B_mass)
{
/*
     This library function assembles the B_mass matrix in
     [B_mass trans][B_mass] for plate elements.

     The specifications for the plate element mass matrix
     was taken from the AMES 232B Winter 1997 notes by
     David Benson, page 79.

		Updated 7/31/06
*/

        int i,i2,i2m1,i2m2;
        for( i = 0; i < npel; ++i )
        {
                i2      =ndof*i+2;
                i2m1    =i2-1;
                i2m2    =i2-2;

                *(B_mass+i2)           = *(shg+i);
                *(B_mass+neqel*1+i2m1) = -*(shg+i);
                *(B_mass+neqel*2+i2m2) = *(shg+i);
        }
        return 1;
}

