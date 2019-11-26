/*
    This utility function assembles either the lumped sum global
    diagonal Mass matrix or, in the case of consistent mass, all
    the element mass matrices.  This is for a finite element program
    which does analysis on a triangle.  It is for modal analysis.

		Updated 11/19/01

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
#include "trconst.h"
#include "trstruct.h"

#define  DEBUG         0

extern int dof, numel, numnp, sof;
extern double shg[sosh], shl[sosh], w[num_int], *Area0;
extern int consistent_mass_flag, consistent_mass_store, lumped_mass_flag;

int matXT(double *, double *, double *, int, int, int);

int triangleB_mass(double *,double *);

int trMassemble(int *connect, double *coord, int *el_matl, int *id,
	double *mass, MATL *matl) 
	
{
        int i, i1, i2, i3, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, node, counter;
	int matl_num;
	double rho, fdum;
        double B_mass[MsoB], B2_mass[MsoB];
        double M_temp[neqlsq], M_el[neqlsq];
        double coord_el_trans[neqel], X1, X2, X3, Y1, Y2, Y3;
        double det[num_int], area_el, wXdet;
        double mass_el[neqel];

/*      initialize all variables  */

        memset(B_mass,0,MsoB*sof);
        memset(B2_mass,0,MsoB*sof);

	memcpy(shg,shl,sosh*sizeof(double));

        for( k = 0; k < numel; ++k )
        {
                matl_num = *(el_matl+k);
                rho = matl[matl_num].rho;
		area_el = 0.0;

/* Zero out the Element mass matrices */

        	memset(M_el,0,neqlsq*sof);
        	memset(mass_el,0,neqel*sof);

/* Create the coord_el transpose vector for one element */

                for( j = 0; j < npel; ++j )
                {
			node = *(connect+npel*k+j);

			*(sdof_el+nsd*j) = nsd*node;
			*(sdof_el+nsd*j+1) = nsd*node+1;

			*(coord_el_trans+j)=*(coord+*(sdof_el+nsd*j));
			*(coord_el_trans+npel*1+j)=*(coord+*(sdof_el+nsd*j+1));

			*(dof_el+ndof*j) = ndof*node;
			*(dof_el+ndof*j+1) = ndof*node+1;
                }

/* Assembly of the Mass matrix. 
*/

		X1 = *(coord_el_trans);
		X2 = *(coord_el_trans + 1);
		X3 = *(coord_el_trans + 2);

		Y1 = *(coord_el_trans + npel*1);
		Y2 = *(coord_el_trans + npel*1 + 1);
		Y3 = *(coord_el_trans + npel*1 + 2);

/*
     It should be noted that I have permutated the node sequences.
     (see Hughes Figure 3.I.5 on page 167; node 3 in Hughes is node 1
     in SLFFEA, node 2 is node 3, and node 1 goes to node 2.)
     This is because I want to be consistant with the tetrahedron.  You can
     read more about this change in teshl for the tetrahedron.

     This change only effects the calculation of the area.  The [M_el] matrix
     is still the same.
*/

		fdum = (X2 - X1)*(Y3 - Y1) - (X3 - X1)*(Y2 - Y1);

		if( fdum <= 0.0 )
		{
			printf("the element (%d) is inverted; 2*Area:%f\n", k, fdum);
			return 0;
		}

/* A factor of 0.5 is needed to do the integration.  See Eq. 3.I.34 in
   "The Finite Element Method" by Thomas Hughes, page 174
*/
		area_el = .5*fdum;

		fdum = area_el*rho/12.0;

#if !DEBUG
/*
   This is taken from "Theory of Matrix Structural Analysis" by
   J. S. Przemieniecki, page 297-299.
*/
		*(M_el)    = 2.0*fdum; *(M_el+2)  = 1.0*fdum; *(M_el+4)  = 1.0*fdum;
		*(M_el+7)  = 2.0*fdum; *(M_el+9)  = 1.0*fdum; *(M_el+11) = 1.0*fdum;
		*(M_el+12) = 1.0*fdum; *(M_el+14) = 2.0*fdum; *(M_el+16) = 1.0*fdum;
		*(M_el+19) = 1.0*fdum; *(M_el+21) = 2.0*fdum; *(M_el+23) = 1.0*fdum;
		*(M_el+24) = 1.0*fdum; *(M_el+26) = 1.0*fdum; *(M_el+28) = 2.0*fdum;
		*(M_el+31) = 1.0*fdum; *(M_el+33) = 1.0*fdum; *(M_el+35) = 2.0*fdum;
#endif

#if DEBUG

/* Assembly of the local mass matrix:

   Below, I am assembling the above mass matrix based on numerical integration.
   The reasons I am doing this are

      1) To help debug the above code
      2) To illustrate 3 point triangle integration.

   Because it is less efficient than simply using a pre-integrated matrix, it has
   been commented out.
*/

		for( j = 0; j < num_int; ++j )
		{

/* Assembly of the B matrix for mass */

		    check = triangleB_mass((shg+npel*(nsd+1)*j + npel*(nsd)),B_mass);
		    if(!check) printf( "Problems with triangleB_mass \n");

/*
		    for( i1 = 0; i1 < nsd; ++i1 )
		    {
			for( i2 = 0; i2 < neqel; ++i2 )
			{
				printf("%16.8e ",*(B_mass+neqel*i1+i2));
			}
			printf(" \n");
		    }
		    printf(" \n");
*/

		    memcpy(B2_mass,B_mass,MsoB*sizeof(double));

		    check=matXT(M_temp, B_mass, B2_mass, neqel, neqel, nsd);
		    if(!check) printf( "Problems with matXT \n");

		    fdum = *(w+j)*rho*area_el;
		    for( i2 = 0; i2 < neqlsq; ++i2 )
		    {
			*(M_el+i2) += *(M_temp+i2)*fdum;
		    }
		}
#endif

/*
		for( i1 = 0; i1 < neqel; ++i1 )
		{
			for( i2 = 0; i2 < neqel; ++i2 )
			{
				printf("%16.8e ",*(M_el+neqel*i1+i2));
			}
			printf(" \n");
		}
		printf(" \n");
*/

                /* printf("This is 2 X Area %10.6f for element %4d\n",2.0*area_el,k);*/
		printf("%4i %16.8e\n",k, area_el);

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
		    printf("This is Area2 %16.8e\n\n",fdum/2.0);

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
