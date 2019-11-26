/*
    This utility function calculates the Area of each
    quad element using shape functions and gaussian
    integration.

                Updated 8/21/06

    SLFFEA source file
    Version:  1.4
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.

*/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qdconst.h"

extern int numel, flag_3D, sof;
extern double shg[sosh], shl[sosh], w[num_int];

int qdshg( double *, int, double *, double *, double *);

int qdArea( int *connect, double *coord, double *Area)
{
	int i, i2, j, k, node, check;
	double fdum1, fdum2, fdum3, fdum4;
	double rotate[npel*nsd2*nsd];
	double coord_el[npel*nsd], coord_el_trans[npel*nsd],
		coord_el_local[npel*nsd2], coord_el_local_trans[npel*nsd2];
	double det[num_int];
	double local_x[nsd], local_y[nsd], local_z[nsd], vec_dum[nsd],
		vec_dum1[nsd], vec_dum2[nsd], vec_dum3[nsd], vec_dum4[nsd],
		xp0[npel], xp1[npel], xp2[npel], yp0[npel], yp1[npel], yp2[npel],
		zp0[npel], zp1[npel], zp2[npel];

	for( k = 0; k < numel; ++k )
	{

/* Create the coord_el transpose vector for one element */

		for( j = 0; j < npel; ++j )
		{
			node=*(connect+npel*k+j);

                        *(coord_el+nsd*j)=*(coord+nsd*node);
                        *(coord_el+nsd*j+1)=*(coord+nsd*node+1);
                        *(coord_el+nsd*j+2)=*(coord+nsd*node+2);

			*(coord_el_trans+j)=*(coord+nsd*node);
			*(coord_el_trans+npel*1+j)=*(coord+nsd*node+1);
			*(coord_el_trans+npel*2+j)=*(coord+nsd*node+2);
		}

		if(!flag_3D)
		{
		    for( j = 0; j < npel; ++j )
		    {
			*(coord_el_local_trans + j) = *(coord_el_trans + j);
			*(coord_el_local_trans + 1*npel + j) = *(coord_el_trans + npel*1 + j);
		    }
		}
		else
		{
/* For 3-D quad meshes, I have to rotate from the global coordinates to the local x and
   y coordinates which lie in the plane of the element.  To do this I have to calculate
   the normal to the plate face, then cross product that normal with an in plane vector.
   Like the shell, I have created a local x and y axis coordinate system at each node.
   I will also be using the algorithm used in my shell code that tries to align the
   normal vector to the global z direction.
*/

		    *(xp2+1) = *(xp1+3) = *(xp0+0) = *(coord_el_trans);
		    *(xp2+2) = *(xp1+0) = *(xp0+1) = *(coord_el_trans + 1);
		    *(xp2+3) = *(xp1+1) = *(xp0+2) = *(coord_el_trans + 2);
		    *(xp2+0) = *(xp1+2) = *(xp0+3) = *(coord_el_trans + 3);

		    *(yp2+1) = *(yp1+3) = *(yp0+0) = *(coord_el_trans + npel*1);
		    *(yp2+2) = *(yp1+0) = *(yp0+1) = *(coord_el_trans + npel*1 + 1);
		    *(yp2+3) = *(yp1+1) = *(yp0+2) = *(coord_el_trans + npel*1 + 2);
		    *(yp2+0) = *(yp1+2) = *(yp0+3) = *(coord_el_trans + npel*1 + 3);

		    *(zp2+1) = *(zp1+3) = *(zp0+0) = *(coord_el_trans + npel*2);
		    *(zp2+2) = *(zp1+0) = *(zp0+1) = *(coord_el_trans + npel*2 + 1);
		    *(zp2+3) = *(zp1+1) = *(zp0+2) = *(coord_el_trans + npel*2 + 2);
		    *(zp2+0) = *(zp1+2) = *(zp0+3) = *(coord_el_trans + npel*2 + 3);

/*
   Calculating rotation matrix for the fiber q[i,j] matrix.
   The algorithm below is taken from "The Finite Element Method" by Thomas Hughes,
   page 388.  The goal is to find the local shell coordinates which come closest
   to the global x, y, z coordinates.  In the algorithm below, vec_dum is set to either
   the global x, y, or z basis vector based on the one of the 2 smaller components of xl.hat
   whose largest component is the local z direction of that node.  Once set, the cross
   product of xl.hat and vec_dum produces the local y fiber direction.  This local y is
   then crossed with xl.hat to produce local x.
*/

		    for( j = 0; j < npel; ++j )
		    {
			*(vec_dum1)     = *(xp1+j) - *(xp0+j);
			*(vec_dum1 + 1) = *(yp1+j) - *(yp0+j);
			*(vec_dum1 + 2) = *(zp1+j) - *(zp0+j);
			*(vec_dum2)     = *(xp2+j) - *(xp0+j);
			*(vec_dum2 + 1) = *(yp2+j) - *(yp0+j);
			*(vec_dum2 + 2) = *(zp2+j) - *(zp0+j);

/* Calculate the local z basis vector for node j */
			check = normcrossX( vec_dum1, vec_dum2, local_z );
			if(!check) printf( "Problems with normcrossX \n");

			fdum1 = fabs(*(local_z));
			fdum2 = fabs(*(local_z+1));
			fdum3 = fabs(*(local_z+2));

			memset(vec_dum,0,nsd*sof);
			i2=1;
			if( fdum1 > fdum3)
			{
			    fdum3=fdum1;
			    i2=2;
			}
			if( fdum2 > fdum3) i2=3;
			*(vec_dum+(i2-1))=1.0;

/* Calculate the local y basis vector for node j */
			check = normcrossX( local_z, vec_dum, local_y );
			if(!check) printf( "Problems with normcrossX \n");

/* Calculate the local x basis vector for node j */
			check = normcrossX( local_y, local_z, local_x );
			if(!check) printf( "Problems with normcrossX \n");

			*(rotate + j*nsd2*nsd) = *(local_x);
			*(rotate + j*nsd2*nsd + 1) = *(local_x + 1);
			*(rotate + j*nsd2*nsd + 2) = *(local_x + 2);
			*(rotate + j*nsd2*nsd + 3) = *(local_y);
			*(rotate + j*nsd2*nsd + 4) = *(local_y + 1);
			*(rotate + j*nsd2*nsd + 5) = *(local_y + 2);

/* Put coord_el into local coordinates */

			check = matX( (coord_el_local+nsd2*j), (rotate + j*nsd2*nsd),
				(coord_el+nsd*j), nsd2, 1, nsd);
			if(!check) printf( "Problems with  matX \n");
			*(coord_el_local_trans + j) = *(coord_el_local+nsd2*j);
			*(coord_el_local_trans + npel*1 + j) = *(coord_el_local+nsd2*j+1);
		    }

		}

/* Assembly of the shg matrix for each integration point */

		check=qdshg(det, k, shl, shg, coord_el_local_trans);
		if(!check) printf( "Problems with qdshg \n");

/* Calculate the Area from determinant of the Jacobian */

		for( j = 0; j < num_int; ++j )
		{
			*(Area + k) += *(w+j)*(*(det+j));
		}
	}

	return 1;
}
