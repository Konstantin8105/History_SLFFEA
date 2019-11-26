/*
     Because the derivatives of triangle elements are constant,
     this file is not really needed.  It is only here for testing
     purposes.

     SLFFEA source file
     Version:  1.3
     Copyright (C) 1999, 2000, 2001, 2002  San Le

     The source code contained in this file is released under the
     terms of the GNU Library General Public License.
*/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "trconst.h"

int dotX(double *,double *, double *, int );

int trshl( int gauss_stress_flag, double *shl, double *w)
{
/* 
     This subroutine calculates the local shape function derivatives for
     a triangle element at the gauss points.

     Because the derivatives of triangle elements are constant,
     this function is not really needed.  It is only here for testing
     purposes.

     It is a C translation of the subroutine QDCSHL from the book
     "The Finite Element Method" by Thomas Hughes, page 784.

     To see the shape functions for triangles, look on page 167.

     Also, look on page 13-3 to 13-6 in the "ANSYS User's Manual".  It
     discusses integration points for triangles.

     It should be noted that I have permutated the node sequences
     given by Hughes so that in terms of the triangle given in
     Figure 3.I.5 on page 167, node 3 in Hughes is node 1 in SLFFEA,
     node 2 is node 3, and node 1 goes to node 2.  This is because
     I want to be consistant with the tetrahedron.  You can read
     more about this change in teshl for the tetrahedron.

 ....  CALCULATE INTEGRATION-RULE WEIGHTS, SHAPE FUNCTIONS AND
       LOCAL DERIVATIVES FOR A EIGHT-NODE HEXAHEDRAL ELEMENT
       r, s = LOCAL ELEMENT COORDINATES along triangle faces
       *(shl+npel*(nsd+1)*k+i) = LOCAL ("r") DERIVATIVE OF SHAPE FUNCTION
       *(shl+npel*(nsd+1)*k+npel*1+i) = LOCAL ("s") DERIVATIVE OF SHAPE FUNCTION
       *(shl+npel*(nsd+1)*k+npel*2+i) = LOCAL SHAPE FUNCTION
       *(w+k)    = INTEGRATION-RULE WEIGHT
       i       = LOCAL NODE NUMBER
       k       = INTEGRATION POINT
       num_int = NUMBER OF INTEGRATION POINTS, EQ. 1 OR 3

			Updated 11/26/09
*/

	double ra[]={ pt1667, pt6667, pt1667};
	double sa[]={ pt1667, pt1667, pt6667};
	double temps,tempr,r,s,fdum;
	int i,j,k;

        if( !gauss_stress_flag )
        {

/* Set ra, sa, ta for shl_node calculation */

                *(ra) = 0.0;     *(sa) = 0.0;
                *(ra + 1) = 1.0; *(sa + 1) = 0.0;
                *(ra + 2) = 0.0; *(sa + 2) = 1.0;
        }

        for( k = 0; k < num_int; ++k )
        {
/* calculating the weights and local dN/ds,dr matrix */

		*(w+k)=pt3333;

        	r=*(ra+k);
        	s=*(sa+k);
		fdum = (1.0 - r - s);

/* dN/dr */
		*(shl+npel*(nsd+1)*k)=   -1.0;         /* node 0 */
		*(shl+npel*(nsd+1)*k+1)=  1.0;         /* node 1 */
		*(shl+npel*(nsd+1)*k+2)=  0.0;	       /* node 2 */

/* dN/ds */
		*(shl+npel*(nsd+1)*k+npel*1)=  -1.0;   /* node 0 */
		*(shl+npel*(nsd+1)*k+npel*1+1)= 0.0;   /* node 1 */
		*(shl+npel*(nsd+1)*k+npel*1+2)= 1.0;   /* node 2 */

/* N */
		*(shl+npel*(nsd+1)*k+npel*2)=   fdum;  /* node 0 */
		*(shl+npel*(nsd+1)*k+npel*2+1)= r;     /* node 1 */
		*(shl+npel*(nsd+1)*k+npel*2+2)= s;     /* node 2 */

	       	/*printf("\n");*/

		/*for( i = 0; i < nsd+1; ++i )
		{
		    for( j = 0; j < npel; ++j )
		    {
			printf(" %14.6e",*(shl+npel*(nsd+1)*k+npel*i + j));
		    }
		    printf("\n");
		}
		printf("\n");*/
	}
        return 1;
}

int trshg( double *det, int el, double *shl, double *shg, double *xl)
{
/*
     This subroutine calculates the global shape function derivatives for
     a quad element at the gauss points.

     Because the derivatives of triangle elements are constant, 
     this function is not really needed.  It is only here for testing
     purposes.

     It is a C translation of the subroutine QDCSHG from the book
     "The Finite Element Method" by Thomas Hughes, page 783.

 ....  CALCULATE GLOBAL DERIVATIVES OF SHAPE FUNCTIONS AND
       JACOBIAN DETERMINANTS FOR A FOUR-NODE QUADRALATERAL ELEMENT

       *(xl+j+npel*i) = GLOBAL COORDINATES(I LOOPS OVER X,Y)
       *(det)  = JACOBIAN DETERMINANT
       *(shl+npel*(nsd+1)*k+i) = LOCAL ("r") DERIVATIVE OF SHAPE FUNCTION
       *(shl+npel*(nsd+1)*k+npel*1+i) = LOCAL ("s") DERIVATIVE OF SHAPE FUNCTION
       *(shl+npel*(nsd+1)*k+npel*2+i) = LOCAL SHAPE FUNCTION
       *(shg+npel*(nsd+1)*k+i) = X-DERIVATIVE OF SHAPE FUNCTION
       *(shg+npel*(nsd+1)*k+npel*1+i) = Y-DERIVATIVE OF SHAPE FUNCTION
       *(shg+npel*(nsd+1)*k+npel*2+i) = shl(npel*(nsd+1)*k+npel*2+i)
       *(xs+2*i+j) = JACOBIAN MATRIX
	  i    = LOCAL NODE NUMBER OR GLOBAL COORDINATE NUMBER
	  j    = GLOBAL COORDINATE NUMBER
	  k    = INTEGRATION-POINT NUMBER
       num_int    = NUMBER OF INTEGRATION POINTS, EQ.3

			Updated 11/20/01
*/

	double xs[4],temp;
	int check,i,j,k;

	memcpy(shg,shl,sosh*sizeof(double));

/* The jacobian dx/dc is calculated below */

#if 0
	for( j = 0; j < nsd; ++j )
	{
		for( i = 0; i < nsd; ++i )
		{
		   check=dotX((xs+nsd*i+j),(shg+npel*j), (xl+npel*i),npel);
		}
	}
#endif

#if 1
	*(xs) =   - *(xl) + *(xl + 1);
	*(xs+1) = - *(xl) + *(xl + 2);
	*(xs+2) = - *(xl + 3) + *(xl + 4);
	*(xs+3) = - *(xl + 3) + *(xl + 5);
#endif

	*(det)=*(xs)*(*(xs+3))-*(xs+2)*(*(xs+1));
	/*printf(" %9.5f %9.5f\n %9.5f %9.5f\n",
		*(xs),*(xs+1),*(xs+nsd*1),*(xs+nsd*1+1));
	printf("%d %f\n", k, *(det));*/

	if(*(det) <= 0.0 )
	{
		printf("the element (%d) is inverted; det:%f\n", el,*(det));
		return 0;
	}

/* The inverse of the jacobian, dc/dx, is calculated below */

	temp=*(xs);
	*(xs)=*(xs+3)/(*(det));
	*(xs+1)*=-1./(*(det));
	*(xs+2)*=-1./(*(det));
	*(xs+3)=temp/(*(det));

	for( k = 0; k < num_int; ++k )
	{
	   for( i = 0; i < npel; ++i )
	   {
		*(shg+npel*(nsd+1)*k+i) = *(xs)*(*(shl+npel*(nsd+1)*k+i))+
			*(xs+2)*(*(shl+npel*(nsd+1)*k+npel*1+i));
		*(shg+npel*(nsd+1)*k+npel*1+i)=*(xs+1)*(*(shl+npel*(nsd+1)*k+i))+
			*(xs+3)*(*(shl+npel*(nsd+1)*k+npel*1+i));
		/*printf("%d %f %f %f\n", i, *(shg+npel*(nsd+1)*k+i),
			*(shg+npel*(nsd+1)*k+npel*1+i),
			*(shg+npel*(nsd+1)*k+npel*2+i));*/
	   }
	}

	return 1;
}


