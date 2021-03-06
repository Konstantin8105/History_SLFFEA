/*
    This utility function assembles either the lumped sum global
    diagonal Mass matrix or, in the case of consistent mass, all
    the element mass matrices.  This is for a finite element program
    which does analysis on a quad.  It is for modal analysis.

		Updated 10/20/06

    SLFFEA source file
    Version:  1.5
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qdconst.h"
#include "qdstruct.h"

extern int dof, numel, numnp, sof, flag_3D, local_vec_flag;
extern double shg[sosh], shl[sosh], w[num_int], *Area0;
extern int consistent_mass_flag, consistent_mass_store, lumped_mass_flag;

int matXrot2(double *, double *, double *, int, int);

int rotXmat2(double *, double *, double *, int, int);

int rotTXmat2(double *, double *, double *, int, int);

int matXrot3(double *, double *, double *, int, int);

int rotXmat3(double *, double *, double *, int, int);

int rotTXmat3(double *, double *, double *, int, int);

int matXT(double *, double *, double *, int, int, int);

int quadB_mass(double *,double *);

int qdshg_mass( double *, int, double *, double *);

int qdMassemble(int *connect, double *coord, int *el_matl, int *id, double *local_xyz,
	double *mass, MATL *matl) 
	
{
	int i, i1, i2, i3, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, node, counter, dum;
	int matl_num;
	double rho, thickness, fdum, fdum1, fdum2, fdum3, fdum4;
	double B_mass[MsoB], B2_mass[MsoB];
	double M_temp[neqlsq], M_el[neqlsq];
	double rotate[npel*nsd2*nsd];
	double coord_el[npel*nsd], coord_el_trans[npel*nsd],
		coord_el_local[npel*nsd2], coord_el_local_trans[npel*nsd2];
	double det[num_int], area_el, wXdet;
	double mass_el[neqel];
	double local_x[nsd], local_y[nsd], local_z[nsd], vec_dum[nsd],
		vec_dum1[nsd], vec_dum2[nsd], vec_dum3[nsd], vec_dum4[nsd];
	double xp[npel], xp0[npel], xp1[npel], xp2[npel],
		yp[npel], yp0[npel], yp1[npel], yp2[npel],
		zp[npel], zp0[npel], zp1[npel], zp2[npel];

/*      initialize all variables  */

	memcpy(shg,shl,sosh*sizeof(double));

	for( k = 0; k < numel; ++k )
	{
		matl_num = *(el_matl+k);
		thickness = matl[matl_num].thick;
		rho = matl[matl_num].rho;
		area_el = 0.0;

/* Create the coord_el transpose vector for one element */

		for( j = 0; j < npel; ++j )
		{
			node = *(connect+npel*k+j);

			*(sdof_el+nsd*j) = nsd*node;
			*(sdof_el+nsd*j+1) = nsd*node+1;
			*(sdof_el+nsd*j+2) = nsd*node+2;

			*(coord_el+nsd*j)=*(coord+*(sdof_el+nsd*j));
			*(coord_el+nsd*j+1)=*(coord+*(sdof_el+nsd*j+1));
			*(coord_el+nsd*j+2)=*(coord+*(sdof_el+nsd*j+2));

			*(coord_el_trans+j)=*(coord+*(sdof_el+nsd*j));
			*(coord_el_trans+npel*1+j)=*(coord+*(sdof_el+nsd*j+1));
			*(coord_el_trans+npel*2+j)=*(coord+*(sdof_el+nsd*j+2));

			*(dof_el+ndof*j) = ndof*node;
			*(dof_el+ndof*j+1) = ndof*node+1;
			*(dof_el+ndof*j+2) = ndof*node+2;
		}

		memset(rotate,0,npel*nsd2*nsd*sof);

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
   y coordinates which lie in the plane of the element.  The local basis used for the
   rotation has already been calculated and stored in local_xyz[].  Below, it is
   copied to rotate[].

   As can be seen below, there are 2 ways to do the rotation.  The first, implemented when:

       local_vec_flag = 0

   has a basis for all four nodes of every element.  The second, implemented when:

       local_vec_flag = 1

   is based on ANSYS where there is only one basis for the element.  You can read more
   about it in:

       ~/slffea-1.4/common/local_vec.c
*/

		    if(!local_vec_flag)
		    {
/* This is the first method. */

			for( j = 0; j < npel; ++j )
			{
			    *(rotate + j*nsd2*nsd) =
				*(local_xyz + npel*nsdsq*k + nsdsq*j);
			    *(rotate + j*nsd2*nsd + 1) =
				*(local_xyz + npel*nsdsq*k + nsdsq*j + 1);
			    *(rotate + j*nsd2*nsd + 2) =
				*(local_xyz + npel*nsdsq*k + nsdsq*j + 2);
			    *(rotate + j*nsd2*nsd + 3) =
				*(local_xyz + npel*nsdsq*k + nsdsq*j + 3);
			    *(rotate + j*nsd2*nsd + 4) =
				*(local_xyz + npel*nsdsq*k + nsdsq*j + 4);
			    *(rotate + j*nsd2*nsd + 5) =
				*(local_xyz + npel*nsdsq*k + nsdsq*j + 5);

/* Put coord_el into local coordinates */

			    check = matX( (coord_el_local+nsd2*j), (rotate + j*nsd2*nsd),
				(coord_el+nsd*j), nsd2, 1, nsd);
			    if(!check) printf( "Problems with  matX \n");
			    *(coord_el_local_trans + j) = *(coord_el_local+nsd2*j);
			    *(coord_el_local_trans + npel*1 + j) = *(coord_el_local+nsd2*j+1);
			}
		    }
		    else
		    {
/* This is the second method. */

			*(rotate)     = *(local_xyz + nsdsq*k);
			*(rotate + 1) = *(local_xyz + nsdsq*k + 1);
			*(rotate + 2) = *(local_xyz + nsdsq*k + 2);
			*(rotate + 3) = *(local_xyz + nsdsq*k + 3);
			*(rotate + 4) = *(local_xyz + nsdsq*k + 4);
			*(rotate + 5) = *(local_xyz + nsdsq*k + 5);

			dum = nsd*npel;
			check = rotXmat2(coord_el_local, rotate, coord_el, 1, dum);
			if(!check) printf( "Problems with  rotXmat2 \n");
			for( j = 0; j < npel; ++j )
			{
			    *(coord_el_local_trans + j) = *(coord_el_local+nsd2*j);
			    *(coord_el_local_trans + npel*1 + j) = *(coord_el_local+nsd2*j+1);
			}
		    }

		}

/* The call to qdshg_mass is only for calculating the determinent */

		check = qdshg_mass(det, k, shg, coord_el_local_trans);
		if(!check) printf( "Problems with qdshg_mass \n");
#if 0
		for( i1 = 0; i1 < num_int; ++i1 )
		{
		    for( i2 = 0; i2 < npel; ++i2 )
		    {
			printf("%10.6f ",*(shl+npel*(nsd2+1)*i1 + npel*(nsd2) + i2));
		    }
		    printf(" \n");
		}
		printf(" \n");
#endif

/* The loop over j below calculates the 4 points of the gaussian integration
   for several quantities */

/* Zero out the Element mass matrices */

		memset(M_el,0,neqlsq*sof);
		memset(mass_el,0,neqel*sof);

		for( j = 0; j < num_int; ++j )
		{
		    memset(B_mass,0,MsoB*sof);
		    memset(B2_mass,0,MsoB*sof);
		    memset(M_temp,0,neqlsq*sof);

/* Assembly of the B matrix for mass */

		    check = quadB_mass((shg+npel*(nsd2+1)*j + npel*(nsd2)),B_mass);
		    if(!check) printf( "Problems with quadB_mass \n");

#if 0
		    for( i1 = 0; i1 < nsd2; ++i1 )
		    {
			for( i2 = 0; i2 < neqel8; ++i2 )
			{
				printf("%9.6f ",*(B_mass+neqel8*i1+i2));
			}
			printf(" \n");
		    }
		    printf(" \n");
#endif

		    memcpy(B2_mass,B_mass,MsoB*sizeof(double));

		    wXdet = *(w+j)*(*(det+j));

/* Calculate the Area from determinant of the Jacobian */

		    area_el += wXdet;

		    check=matXT(M_temp, B_mass, B2_mass, neqel8, neqel8, nsd2);
		    if(!check) printf( "Problems with matXT \n");

/* Normally, with plane elements, we assume a unit thickness in the transverse direction.  But
   because these elements can be 3 dimensional, I multiply the material property matrix by the
   thickness.  This is justified by equation 5.141a in "Theory of Matrix Structural
   Analysis" by J. S. Przemieniecki on page 86.
*/

		    fdum = rho*thickness*wXdet;
		    for( i2 = 0; i2 < neqlsq; ++i2 )
		    {
			*(M_el+i2) += *(M_temp+i2)*fdum;
		    }
		}

		if(!flag_3D)
		{
/* For 2-D meshes */
		   for( i = 0; i < npel; ++i )
		   {
			for( j = 0; j < npel; ++j )
			{
/* row for displacement x */
			   *(M_temp + ndof*neqel*i + ndof*j) =
				*(M_el + ndof2*neqel8*i + ndof2*j);
			   *(M_temp + ndof*neqel*i + ndof*j + 1) =
				*(M_el + ndof2*neqel8*i + ndof2*j + 1);
			   *(M_temp + ndof*neqel*i + ndof*j + 2) = 0.0;

/* row for displacement y */
			   *(M_temp + ndof*neqel*i + neqel  + ndof*j) =
				*(M_el + ndof2*neqel8*i + neqel8 + ndof2*j);
			   *(M_temp + ndof*neqel*i + neqel + ndof*j + 1) =
				*(M_el + ndof2*neqel8*i + neqel8 + ndof2*j + 1);
			   *(M_temp + ndof*neqel*i + neqel + ndof*j + 2) = 0.0;

/* row for displacement z */
			   *(M_temp + ndof*neqel*i + 2*neqel + ndof*j) = 0.0;
			   *(M_temp + ndof*neqel*i + 2*neqel + ndof*j + 1) = 0.0;
			   *(M_temp + ndof*neqel*i + 2*neqel + ndof*j + 2) = 0.0;
			}
		   }
		   memcpy(M_el, M_temp, neqlsq*sizeof(double));
		}
		else
		{
/* For 3-D meshes */

/* Put M back to global coordinates */

		   if(!local_vec_flag)
		   {
			check = matXrot3(M_temp, M_el, rotate, neqel8, neqel);
			if(!check) printf( "Problems with matXrot3 \n");

			check = rotTXmat3(M_el, rotate, M_temp, neqel, neqel);
			if(!check) printf( "Problems with rotTXmat3 \n");
		   }
		   else
		   {
			check = matXrot2(M_temp, M_el, rotate, neqel8, neqel);
			if(!check) printf( "Problems with matXrot2 \n");

			check = rotTXmat2(M_el, rotate, M_temp, neqel, neqel);
			if(!check) printf( "Problems with rotTXmat2 \n");
		   }
		}

		/* printf("This is 2 X Area %10.6f for element %4d\n",2.0*area_el,k);*/

		if(lumped_mass_flag)
		{

/* Creating the diagonal lumped mass Matrix */

		    fdum = 0.0;
		    for( i2 = 0; i2 < neqel; ++i2 )
		    {   
			/*printf("This is mass_el for el %3d",k);*/
			for( i3 = 0; i3 < neqel; ++i3 )
			{
			    *(mass_el+i2) += *(M_el+neqel*i2+i3);
			}
			/*printf("%9.6f\n\n",*(mass_el+i2));*/
			fdum += *(mass_el+i2);
		    }   
		    /*printf("This is Area2 %9.6f\n\n",fdum);*/

		    for( j = 0; j < neqel; ++j )
		    {   
			*(mass+*(dof_el+j)) += *(mass_el + j);
		    }
		}

		if(consistent_mass_flag)
		{

/* Storing all the element mass matrices */

		    for( j = 0; j < neqlsq; ++j )
		    {   
			*(mass + neqlsq*k + j) = *(M_el + j);
		    }
		}
	}

	if(lumped_mass_flag)
	{
/* Contract the global mass matrix using the id array only if lumped
   mass is used. */

	    counter = 0;
	    for( i = 0; i < dof ; ++i )
	    {
		/* printf("%5d  %16.8e\n", i, *(mass+i));*/
		if( *(id + i ) > -1 )
		{
		    *(mass + counter ) = *(mass + i );
		    ++counter;
		}
	    }
	}

	return 1;
}
