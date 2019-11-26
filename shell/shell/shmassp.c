/*
    This utility function takes the product of a vector with the
    consistent mass matrix.  This is for a finite element program
    which does analysis on a shell.  It is for modal analysis.

		Updated 10/10/06

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

int shMassPassemble(int *connect, double *coord, int *el_matl, double *lamina_ref,
	double *fiber_xyz, double *mass, MATL *matl, double *P_global, double *U) 
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
	double det[num_int];
	XL xl;
	double P_el[neqel];

	memset(P_global,0,dof*sof);

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

			xl.bar[j] = *(lamina_ref + nsd*node);
			xl.bar[npell*1+j] = *(lamina_ref + nsd*node + 1);
			xl.bar[npell*2+j] = *(lamina_ref + nsd*node + 2);

			fdum1=*(coord_el_trans+npell+j)-*(coord_el_trans+j);
			fdum2=*(coord_el_trans+npel*1+npell+j)-*(coord_el_trans+npel*1+j);
			fdum3=*(coord_el_trans+npel*2+npell+j)-*(coord_el_trans+npel*2+j);
			fdum4=sqrt(fdum1*fdum1+fdum2*fdum2+fdum3*fdum3);

			*(zp1+j)=.5*(1.0-zeta)*fdum4;
			*(zm1+j)=-.5*(1.0+zeta)*fdum4;

			xl.hat[j] = *(fiber_xyz + nsdsq*node + 2*nsd);
			xl.hat[npell*1+j] = *(fiber_xyz + nsdsq*node + 2*nsd + 1);
			xl.hat[npell*2+j] = *(fiber_xyz + nsdsq*node + 2*nsd + 2);


/* Create the rotation matrix */

			for( i1 = 0; i1 < nsd; ++i1 )
			{
			    rotate.f_shear[nsdsq*j + i1] =
				*(fiber_xyz + nsdsq*node + i1);
			    rotate.f_shear[nsdsq*j + 1*nsd + i1] =
				*(fiber_xyz + nsdsq*node + 1*nsd + i1);
			    rotate.f_shear[nsdsq*j + 2*nsd + i1] =
				*(fiber_xyz + nsdsq*node + 2*nsd + i1);
			}
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

/* Loop through the integration points */

/* Zero out the Element mass matrices */

		memset(M_el,0,neqlsq*sof);

/* The loop over i4 below calculates the 2 fiber points of the gaussian integration */
		for( i4 = 0; i4 < num_ints; ++i4 )
		{

/* The loop over j below calculates the 2X2 points of the gaussian integration
   over the lamina for several quantities */

		   for( j = 0; j < num_intb; ++j )
		   {
			memset(B_mass,0,MsoB*sof);
			memset(B2_mass,0,MsoB*sof);
			memset(M_temp,0,neqlsq*sof);

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
