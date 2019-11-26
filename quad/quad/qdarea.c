/*
    This utility function calculates the Area of each
    quad element using shape functions and gaussian
    integration.

                Updated 9/27/01

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
#include "qdconst.h"

extern int numel;
extern double shg[sosh], shl[sosh], w[num_int];

int qdshg( double *, int, double *, double *, double *);

int qdArea( int *connect, double *coord, double *Area)
{
	double coord_el_trans[npel*nsd];
	int i, j, k, node, check;
	double det[num_int];

	for( k = 0; k < numel; ++k )
	{

/* Create the coord_el transpose vector for one element */

		for( j = 0; j < npel; ++j )
		{
			node=*(connect+npel*k+j);

			*(coord_el_trans+j)=*(coord+nsd*node);
			*(coord_el_trans+npel*1+j)=*(coord+nsd*node+1);
		}

		check = qdshg(det, k, shl, shg, coord_el_trans);
		if(!check) printf( "Problems with qdshg \n");

/* Calculate the Area from determinant of the Jacobian */

		for( j = 0; j < num_int; ++j )
		{
			*(Area + k) += *(w+j)*(*(det+j));
		}
	}

	return 1;
}
