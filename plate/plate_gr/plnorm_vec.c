/*
    This program Calculates the normal vectors of a mesh
    for plate elements.
  
  		Last Update 3/10/01

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
#include "../plate/plconst.h"
#include "plstrcgr.h"

extern int numel;

int normcrossX(double *, double *, double *);

int plnormal_vectors(int *connecter, double *coord, NORM *norm, double *zcoord)
{
        int i, i2, j, k, sdof_el[npel*nsd], ii, check, counter, node;
	int l,m,n;
        double coord_el[npel*3];
	double d1[3], d2[3], norm_temp[3];

        for( k = 0; k < numel; ++k )
        {
                for( j = 0; j < npel; ++j )
                {

/* Calculate element degrees of freedom */

                        node = *(connecter+npel*k+j);
                        *(sdof_el+nsd*j) = nsd*node;
                        *(sdof_el+nsd*j+1) = nsd*node+1;

/* Calculate local coordinates */

                        *(coord_el+3*j)=*(coord+*(sdof_el+nsd*j));
                        *(coord_el+3*j+1)=*(coord+*(sdof_el+nsd*j+1));
                        *(coord_el+3*j+2)=*(zcoord+node);

    			/*printf( "coord %9.5f %9.5f %9.5f \n",*(coord_el+3*j),
				*(coord_el+3*j+1),*(coord_el+3*j+2));*/
                }

/* Calculate normal vectors
   Note that the normal of the faces are calculated differently from the brick
   and the shell.
*/

	     	*(d1)=*(coord_el+3)-*(coord_el);
	     	*(d1+1)=*(coord_el+4)-*(coord_el+1);
	     	*(d1+2)=*(coord_el+5)-*(coord_el+2);
	     	*(d2)=*(coord_el+9)-*(coord_el);
	     	*(d2+1)=*(coord_el+10)-*(coord_el+1);
	     	*(d2+2)=*(coord_el+11)-*(coord_el+2);
		normcrossX(d1, d2, norm_temp);
		norm[k].face[0].x = *(norm_temp);
		norm[k].face[0].y = *(norm_temp+1);
		norm[k].face[0].z = *(norm_temp+2);
    		/*printf( "%9.5f %9.5f %9.5f \n",*(norm_temp),
			*(norm_temp+1),*(norm_temp+2));*/

	     	*(d1)=*(coord_el+9)-*(coord_el+6);
	     	*(d1+1)=*(coord_el+10)-*(coord_el+7);
	     	*(d1+2)=*(coord_el+11)-*(coord_el+8);
	     	*(d2)=*(coord_el+3)-*(coord_el+6);
	     	*(d2+1)=*(coord_el+4)-*(coord_el+7);
	     	*(d2+2)=*(coord_el+5)-*(coord_el+8);
		normcrossX(d1, d2, norm_temp);
		norm[k].face[1].x = *(norm_temp);
		norm[k].face[1].y = *(norm_temp+1);
		norm[k].face[1].z = *(norm_temp+2);
    		/*printf( "%9.5f %9.5f %9.5f \n",*(norm_temp),
			*(norm_temp+1),*(norm_temp+2));*/
	}
	return 1;
}
