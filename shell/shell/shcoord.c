/*
    This utility function uses the the element thicknesses of a singly curved 4
    node shell element and the normal vectors at each node to calculate the 
    corresponding top nodal coordinates.  It then does an averaging based on the
    number of elements which share a particular node.

		Updated 8/10/06

    SLFFEA source file
    Version:  1.4
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shconst.h"
#include "shstruct.h"

int normcrossX(double *, double *, double *);

int matXT(double *, double *, double *, int, int, int);

extern sdof, numel, numnp;

int shTopCoordinates(int *connect, double *coord, int *el_matl, double *fiber_vec,
	MATL *matl, double *node_counter)
{
	int i, i1, i2, j, k, sdof_el[npel*nsd];
	int check, counter, node;
	int matl_num;
	double fdum1, fdum2, fdum3, fdum4;
	double coord_el_trans[npel*nsd];
	double thickness;

	for( k = 0; k < numel; ++k )
	{

		matl_num = *(el_matl+k);
		thickness = matl[matl_num].thick;

/* Calculate the corresponding top coordinates by first setting it equal to the components
   of the normal node vector multiplied by the element thickness.  Sum this over every element
   which shares the node.
*/

		for( j = 0; j < npell; ++j )
		{
			node = *(connect+npell*k+j);

			*(sdof_el+nsd*j)=nsd*node;
			*(sdof_el+nsd*j+1)=nsd*node+1;
			*(sdof_el+nsd*j+2)=nsd*node+2;

			*(coord+*(sdof_el+nsd*j)+sdof) +=
				*(fiber_vec+*(sdof_el+nsd*j))*thickness;
			*(coord+*(sdof_el+nsd*j+1)+sdof) +=
				*(fiber_vec+*(sdof_el+nsd*j+1))*thickness;
			*(coord+*(sdof_el+nsd*j+2)+sdof) +=
				*(fiber_vec+*(sdof_el+nsd*j+2))*thickness;
		}

	}

/* Average all the thickness components of the top nodes and then add it to the
   bottom nodes */

	for( i = 0; i < numnp ; ++i )
	{
		*(coord+nsd*i+sdof) /= *(node_counter+i);
		*(coord+nsd*i+1+sdof) /= *(node_counter+i);
		*(coord+nsd*i+2+sdof) /= *(node_counter+i);

		*(coord+nsd*i+sdof) += *(coord+nsd*i);
		*(coord+nsd*i+1+sdof) += *(coord+nsd*i+1);
		*(coord+nsd*i+2+sdof) += *(coord+nsd*i+2);

		/*printf("\n node %3d %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
                        i,*(coord+nsd*i+sdof),*(coord+nsd*i+1+sdof),*(coord+nsd*i+2+sdof),
                        *(coord+nsd*i),*(coord+nsd*i+1),*(coord+nsd*i+2));*/

	}

	return 1;
}

