/*
    This utility function calculates the normal vectors for each node on a
    suface mesh made up of 4 node elements.  This is done element by element
    where the normal is set to the cross product of the 2 adjacent sides which
    share a node.  Because there may be multiple elements which may contain
    the node, averaging is done at the end based on the number of elements which
    have the node.

		Updated 5/2/05

    SLFFEA source file
    Version:  1.4
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define nsd        3            /* spatial dimensions per node */
#define npel       4            /* nodes per element */

int normal(double *);

int normcrossX(double *, double *, double *);

int matXT(double *, double *, double *, int, int, int);

int node_normals(int *connect, double *coord, double *norm_vec, double *node_counter,
	int numel, int numnp)
{
	int i, i1, i2, j, k, sdof_el[npel*nsd];
	int check, counter, node;
	int matl_num;
	double fdum1, fdum2, fdum3, fdum4;
	double coord_el_trans[npel*nsd];
	double vec_dum1[nsd], vec_dum2[nsd], vec_dum3[nsd];

	for( k = 0; k < numel; ++k )
	{

/* Create the coord transpose vector and other variables for one element */

		for( j = 0; j < npel; ++j )
		{
			node = *(connect+npel*k+j);

			*(sdof_el+nsd*j)=nsd*node;
			*(sdof_el+nsd*j+1)=nsd*node+1;
			*(sdof_el+nsd*j+2)=nsd*node+2;

			*(coord_el_trans+j) =
				*(coord+*(sdof_el+nsd*j));
			*(coord_el_trans+npel*1+j) =
				*(coord+*(sdof_el+nsd*j+1));
			*(coord_el_trans+npel*2+j) =
				*(coord+*(sdof_el+nsd*j+2));

/* Count the number of times a particular node is part of an element */

			*(node_counter+node) += 1.0;

		}

/* Calculate the normal node vectors */

/* for node 0 */
		*(vec_dum1) = *(coord_el_trans+1) - *(coord_el_trans);
		*(vec_dum1+1) = *(coord_el_trans+npel+1) - *(coord_el_trans+npel);
		*(vec_dum1+2) = *(coord_el_trans+2*npel+1) - *(coord_el_trans+2*npel);

		*(vec_dum2) = *(coord_el_trans+3) - *(coord_el_trans);
		*(vec_dum2+1) = *(coord_el_trans+npel+3) - *(coord_el_trans+npel);
		*(vec_dum2+2) = *(coord_el_trans+2*npel+3) - *(coord_el_trans+2*npel);

		/*printf("vvvv %14.5e %14.5e %14.5e\n", *(vec_dum1), *(vec_dum1+1), *(vec_dum1+2));
		printf("vvvv %14.5e %14.5e %14.5e\n", *(vec_dum2), *(vec_dum2+1), *(vec_dum2+2));*/
		check = normcrossX(vec_dum1, vec_dum2, vec_dum3);
		if(!check) printf( "Problems with normcrossX \n");
		*(norm_vec+*(sdof_el)) += *(vec_dum3);
		*(norm_vec+*(sdof_el+1)) += *(vec_dum3+1);
		*(norm_vec+*(sdof_el+2)) += *(vec_dum3+2);
		/*printf("1111 %14.5e %14.5e %14.5e\n", *(norm_vec+*(sdof_el)),
			*(norm_vec+*(sdof_el+1)),
			*(norm_vec+*(sdof_el+2)));*/
/* for node 1 */
		*(vec_dum1) = *(coord_el_trans+2) - *(coord_el_trans+1);
		*(vec_dum1+1) = *(coord_el_trans+npel+2) - *(coord_el_trans+npel+1);
		*(vec_dum1+2) = *(coord_el_trans+2*npel+2) - *(coord_el_trans+2*npel+1);

		*(vec_dum2) = *(coord_el_trans) - *(coord_el_trans+1);
		*(vec_dum2+1) = *(coord_el_trans+npel) - *(coord_el_trans+npel+1);
		*(vec_dum2+2) = *(coord_el_trans+2*npel) - *(coord_el_trans+2*npel+1);

		check = normcrossX(vec_dum1, vec_dum2, vec_dum3);
		if(!check) printf( "Problems with normcrossX \n");
		*(norm_vec+*(sdof_el+nsd*1)) += *(vec_dum3);
		*(norm_vec+*(sdof_el+nsd*1+1)) += *(vec_dum3+1);
		*(norm_vec+*(sdof_el+nsd*1+2)) += *(vec_dum3+2);
/* for node 2 */
		*(vec_dum1) = *(coord_el_trans+3) - *(coord_el_trans+2);
		*(vec_dum1+1) = *(coord_el_trans+npel+3) - *(coord_el_trans+npel+2);
		*(vec_dum1+2) = *(coord_el_trans+2*npel+3) - *(coord_el_trans+2*npel+2);

		*(vec_dum2) = *(coord_el_trans+1) - *(coord_el_trans+2);
		*(vec_dum2+1) = *(coord_el_trans+npel+1) - *(coord_el_trans+npel+2);
		*(vec_dum2+2) = *(coord_el_trans+2*npel+1) - *(coord_el_trans+2*npel+2);

		check = normcrossX(vec_dum1, vec_dum2, vec_dum3);
		if(!check) printf( "Problems with normcrossX \n");
		*(norm_vec+*(sdof_el+nsd*2)) += *(vec_dum3);
		*(norm_vec+*(sdof_el+nsd*2+1)) += *(vec_dum3+1);
		*(norm_vec+*(sdof_el+nsd*2+2)) += *(vec_dum3+2);
/* for node 3 */
		*(vec_dum1) = *(coord_el_trans) - *(coord_el_trans+3);
		*(vec_dum1+1) = *(coord_el_trans+npel) - *(coord_el_trans+npel+3);
		*(vec_dum1+2) = *(coord_el_trans+2*npel) - *(coord_el_trans+2*npel+3);

		*(vec_dum2) = *(coord_el_trans+2) - *(coord_el_trans+3);
		*(vec_dum2+1) = *(coord_el_trans+npel+2) - *(coord_el_trans+npel+3);
		*(vec_dum2+2) = *(coord_el_trans+2*npel+2) - *(coord_el_trans+2*npel+3);

		check = normcrossX(vec_dum1, vec_dum2, vec_dum3);
		if(!check) printf( "Problems with normcrossX \n");
		*(norm_vec+*(sdof_el+nsd*3)) += *(vec_dum3);
		*(norm_vec+*(sdof_el+nsd*3+1)) += *(vec_dum3+1);
		*(norm_vec+*(sdof_el+nsd*3+2)) += *(vec_dum3+2);

	}

/* Average all the components of norm_vec at the nodes */

	for( i = 0; i < numnp ; ++i )
	{
		*(norm_vec+nsd*i) /= *(node_counter+i);
		*(norm_vec+nsd*i+1) /= *(node_counter+i);
		*(norm_vec+nsd*i+2) /= *(node_counter+i);
		check = normal((norm_vec+nsd*i));
		if(!check) printf( "Problems with normal \n");
	}

	return 1;
}

