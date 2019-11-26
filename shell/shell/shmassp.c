/*
    This utility function takes the product of a vector with the
    consistent mass matrix.  This is for a finite element program
    which does analysis on a shell.  It is for modal analysis.

		Updated 9/26/01

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001  San Le 

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

extern int numel, numnp, dof, sof;
extern SH shg, shl;
extern ROTATE rotate, rotate_node;
extern double w[num_int], *Vol0;
extern int consistent_mass_flag, consistent_mass_store;

int matXT(double *, double *, double *, int, int, int);

int matX(double *, double *, double *, int, int, int);

int shellB_mass(double *, double *, double *, double *);

int shshg_mass( double *, int , SH , XL , double *, double *, double *,
        double *);

int normcrossX(double *, double *, double *);

int shMassPassemble(int *connect, double *coord, int *el_matl, double *mass,
	MATL *matl, double *P_global, double *U) 
{
        int i, i1, i2, i3, i4, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, node, counter;
	int matl_num;
	double rho, fdum, fdum1, fdum2, fdum3, fdum4;
        double B_mass[MsoB], B2_mass[MsoB];
        double M_temp[neqlsq], M_el[neqlsq];
	double U_el[neqel];
	double coord_el_trans[npel*nsd], zm1[npell], zp1[npell],
		znode[npell*num_ints], dzdt_node[npell];
	double vec_dum[nsd];
        double det[num_int];
	XL xl;
	double P_el[neqel];

	memset(P_global,0,dof*sof);
	memset(B_mass,0,MsoB*sof);
	memset(B2_mass,0,MsoB*sof);

	if(consistent_mass_store)
	{

/* Assemble P matrix using stored element mass matrices */

	    for( k = 0; k < numel; ++k )
	    {

		for( j = 0; j < npell; ++j )
		{
			node = *(connect+npell*k+j);

			*(dof_el+ndof*j) = ndof*node;
			*(dof_el+ndof*j+1) = ndof*node+1;
			*(dof_el+ndof*j+2) = ndof*node+2;
			*(dof_el+ndof*j+3) = ndof*node+3;
			*(dof_el+ndof*j+4) = ndof*node+4;
		}

/* Assembly of the global P matrix */

		for( j = 0; j < neqel; ++j )
		{
			*(U_el + j) = *(U + *(dof_el+j));
		}

		check = matX(P_el, (mass + k*neqlsq), U_el, neqel, 1, neqel);
		if(!check) printf( "Problems with matX \n");

		for( j = 0; j < neqel; ++j )
		{
			*(P_global+*(dof_el+j)) += *(P_el+j);
		}
	    }
	}
	else
	{

/* Assemble P matrix by re-deriving element mass matrices */

            for( k = 0; k < numel; ++k )
            {
                matl_num = *(el_matl+k);
                rho = matl[matl_num].rho;

/* Zero out the Element mass matrices */

        	memset(M_el,0,neqlsq*sof);

/* Create the coord_el transpose vector and other variables for one element */

                for( j = 0; j < npell; ++j )
                {
			node = *(connect+npell*k+j);

                        *(sdof_el+nsd*j)=nsd*node;
                        *(sdof_el+nsd*j+1)=nsd*node+1;
                        *(sdof_el+nsd*j+2)=nsd*node+2;

                        *(sdof_el+nsd*npell+nsd*j)=nsd*(node+numnp);
                        *(sdof_el+nsd*npell+nsd*j+1)=nsd*(node+numnp)+1;
                        *(sdof_el+nsd*npell+nsd*j+2)=nsd*(node+numnp)+2;

                        *(dof_el+ndof*j) = ndof*node;
                        *(dof_el+ndof*j+1) = ndof*node+1;
                        *(dof_el+ndof*j+2) = ndof*node+2;
                        *(dof_el+ndof*j+3) = ndof*node+3;
                        *(dof_el+ndof*j+4) = ndof*node+4;

/* Create the coord -/+*/

                        *(coord_el_trans+j) =
				*(coord+*(sdof_el+nsd*j));
                        *(coord_el_trans+npel*1+j) =
				*(coord+*(sdof_el+nsd*j+1));
                        *(coord_el_trans+npel*2+j) =
				*(coord+*(sdof_el+nsd*j+2));

                        *(coord_el_trans+npell+j) =
				*(coord+*(sdof_el+nsd*npell+nsd*j));
                        *(coord_el_trans+npel*1+npell+j) =
				*(coord+*(sdof_el+nsd*npell+nsd*j+1));
                        *(coord_el_trans+npel*2+npell+j) =
				*(coord+*(sdof_el+nsd*npell+nsd*j+2));

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
   Calculating rotation matrix for the fiber q[i,j] matrix.

   The algorithm below is taken from "The Finite Element Method" by Thomas Hughes,
   page 388.  The goal is to find the local shell coordinates which come closest
   to the global x, y, z coordinates.  In the algorithm below, vec_dum is set to either
   the global x, y, or z basis vector based on the one of the 2 smaller components of
   xl.hat.  xl.hat itself is the local z direction of that node.  Once set, the cross
   product of xl.hat and vec_dum produces the local y fiber direction.  This local y is
   then crossed with xl.hat to produce local x.
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

/* Set the local z basis vector = xl.hat */
                	rotate.f_shear[nsdsq*j+2*nsd]=xl.hat[j];
                	rotate.f_shear[nsdsq*j+2*nsd+1]=xl.hat[npell*1+j];
                	rotate.f_shear[nsdsq*j+2*nsd+2]=xl.hat[npell*2+j];

/* Calculate the local y basis vector = local z X vec_dum */
                	check = normcrossX((rotate.f_shear+nsdsq*j+2*nsd),
			    vec_dum,(rotate.f_shear+nsdsq*j+1*nsd));
                	if(!check) printf( "Problems with normcrossX \n");

/* Calculate the local x basis vector = local y X local z */
                	check = normcrossX((rotate.f_shear+nsdsq*j+1*nsd),
			    (rotate.f_shear+nsdsq*j+2*nsd),(rotate.f_shear+nsdsq*j));
                	if(!check) printf( "Problems with normcrossX \n");

                }

/* The call to shshg_mass is only for calculating the determinent */

		check = shshg_mass( det, k, shl, xl, zp1, zm1, znode,
			dzdt_node);
		if(!check) printf( "Problems with shshg_mass \n");

#if 0
                for( i1 = 0; i1 < num_intb; ++i1 )
                {
                    for( i2 = 0; i2 < npell; ++i2 )
                    {
                    	printf("%10.6f ",*(shl+npell*(nsd+1)*i1 + npell*(nsd) + i2));
                    }
                    printf(" \n");
                }
                printf(" \n");
#endif

/* The loop over i4 below calculates the 2 fiber points of the gaussian integration */
		for( i4 = 0; i4 < num_ints; ++i4 )
		{

/* The loop over j below calculates the 2X2 points of the gaussian integration
   over the lamina for several quantities */

		   for( j = 0; j < num_intb; ++j )
		   {

/* Assembly of the B matrix for mass */

			check =
			    shellB_mass((shg.bend+npell*(nsd+1)*num_intb*i4+npell*(nsd+1)*j+npell*(nsd)),
			    znode, B_mass, rotate.f_shear);
			if(!check) printf( "Problems with shellB_mass \n");
#if 0
                   	for( i1 = 0; i1 < nsd; ++i1 )
                   	{
                    	    for( i2 = 0; i2 < neqel; ++i2 )
                    	    {
                        	printf("%9.6f ",*(B_mass+neqel*i1+i2));
                    	    }
                    	    printf(" \n");
                   	}
                   	printf(" \n");
#endif

			memcpy(B2_mass,B_mass,MsoB*sizeof(double));

                	check=matXT(M_temp, B_mass, B2_mass, neqel, neqel, nsd);
                	if(!check) printf( "Problems with matXT \n");

			fdum =  rho*(*(w+num_intb*i4+j))*(*(det+num_intb*i4+j));
			for( i2 = 0; i2 < neqlsq; ++i2 )
			{
				*(M_el+i2) += *(M_temp+i2)*fdum;
			}
		   }
		}

/* Assembly of the global P matrix */

		for( j = 0; j < neqel; ++j )
		{
			*(U_el + j) = *(U + *(dof_el+j));
		}

		check = matX(P_el, M_el, U_el, neqel, 1, neqel);
		if(!check) printf( "Problems with matX \n");

		for( j = 0; j < neqel; ++j )
		{
			*(P_global+*(dof_el+j)) += *(P_el+j);
		}

	    }
	}

        return 1;
}
