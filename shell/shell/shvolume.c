/*
    This utility function calculates the Volume of each
    shell element using shape functions and gaussian
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
#include "shconst.h"
#include "shstruct.h"

extern int numel, numnp, sof;
extern SH shg, shl;
extern ROTATE rotate;
extern double w[num_int];

int shshg( double *, int , SH , SH , XL , double *, double *, double *,
	double *, ROTATE );

int shVolume( int *connect, double *coord, double *Vol)
{
        int i, i1, i2, i3, i4, i5, j, k;
	int check, node;
	int matl_num;
	double fdum1, fdum2, fdum3, fdum4;
        double coord_el_trans[npel*nsd], zm1[npell], zp1[npell],
		znode[npell*num_ints], dzdt_node[npell];
        double vec_dum[nsd];
        double det[num_int+num_ints];
	XL xl;

        for( k = 0; k < numel; ++k )
        {

/* Create the coord transpose vector and other variables for one element */

                for( j = 0; j < npell; ++j )
                {
			node = *(connect+npell*k+j);

/* Create the coord -/+*/

                        *(coord_el_trans+j) = *(coord+nsd*node);
                        *(coord_el_trans+npel*1+j) = *(coord+nsd*node+1);
                        *(coord_el_trans+npel*2+j) = *(coord+nsd*node+2);

                        *(coord_el_trans+npell+j) = *(coord+nsd*(node+numnp));
                        *(coord_el_trans+npel*1+npell+j) = *(coord+nsd*(node+numnp)+1);
                        *(coord_el_trans+npel*2+npell+j) = *(coord+nsd*(node+numnp)+2);

/* Create the coord_bar and coord_hat vector for one element */

                        xl.bar[j]=.5*( *(coord_el_trans+j)*(1.0-zeta)+
				*(coord_el_trans+npell+j)*(1.0+zeta));
                        xl.bar[npell*1+j]=.5*( *(coord_el_trans+npel*1+j)*(1.0-zeta)+
				*(coord_el_trans+npel*1+npell+j)*(1.0+zeta));
                        xl.bar[npell*2+j]=.5*( *(coord_el_trans+npel*2+j)*(1.0-zeta)+
				*(coord_el_trans+npel*2+npell+j)*(1.0+zeta));

                        xl.hat[j]=*(coord_el_trans+npell+j)-*(coord_el_trans+j);
                        xl.hat[npell*1+j]=*(coord_el_trans+npel*1+npell+j)-
				*(coord_el_trans+npel*1+j);
                        xl.hat[npell*2+j]=*(coord_el_trans+npel*2+npell+j)-
				*(coord_el_trans+npel*2+j);

			fdum1=fabs(xl.hat[j]);
			fdum2=fabs(xl.hat[npell*1+j]);
			fdum3=fabs(xl.hat[npell*2+j]);
			fdum4=sqrt(fdum1*fdum1+fdum2*fdum2+fdum3*fdum3);
                        xl.hat[j] /= fdum4;
                        xl.hat[npell*1+j] /= fdum4;
                        xl.hat[npell*2+j] /= fdum4;
			*(zp1+j)=.5*(1.0-zeta)*fdum4;
			*(zm1+j)=-.5*(1.0+zeta)*fdum4;
/*
   calculating rotation matrix for the fiber q[i,j] matrix 
*/
        	       	memset(vec_dum,0,nsd*sof);
			i2=1;
			if( fdum1 > fdum3)
			{
				fdum3=fdum1;
				i2=2;
			}
			if( fdum2 > fdum3) i2=3;
			*(vec_dum+(i2-1))=1.0;
                	rotate.f_shear[nsdsq*j+2*nsd]=xl.hat[j];
                	rotate.f_shear[nsdsq*j+2*nsd+1]=xl.hat[npell*1+j];
                	rotate.f_shear[nsdsq*j+2*nsd+2]=xl.hat[npell*2+j];
                	check = normcrossX((rotate.f_shear+nsdsq*j+2*nsd),
			    vec_dum,(rotate.f_shear+nsdsq*j+1*nsd));
                	if(!check) printf( "Problems with normcrossX \n");
                	check = normcrossX((rotate.f_shear+nsdsq*j+1*nsd),
			    (rotate.f_shear+nsdsq*j+2*nsd),(rotate.f_shear+nsdsq*j));
                	if(!check) printf( "Problems with normcrossX \n");
                }

/* Assembly of the shg matrix for each integration point */

		check=shshg( det, k, shl, shg, xl, zp1, zm1, znode,
			dzdt_node, rotate);
		if(!check) printf( "Problems with shshg \n");

/* Calculate the Volume from determinant of the Jacobian */

		for( i4 = 0; i4 < num_ints; ++i4 )
                {
                   for( j = 0; j < num_intb; ++j )
                   {
			*(Vol + k) += *(w+num_intb*i4+j)*(*(det+num_intb*i4+j));

                   }
                }

	}

	return 1;
}

