/*
    This utility function assembles either the lumped sum global
    diagonal Mass matrix or, in the case of consistent mass, all
    the element mass matrices.  This is for a finite element program
    which does analysis on a plate.  It is for modal analysis.

		Updated 9/15/00

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "plconst.h"
#include "plstruct.h"

extern int dof, numel, numnp, sof;
extern SH shg, shl;
extern double w[num_int+1], *Area0;
extern int consistent_mass_flag, consistent_mass_store, lumped_mass_flag;

int matXT(double *, double *, double *, int, int, int);

int plateB_mass(double *,double *);

int plshg_mass( double *, int, SH, double *, double *);

int plMassemble(int *connect, double *coord, int *el_matl, int *id,
	double *mass, MATL *matl) 
	
{
        int i, i1, i2, i3, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, node, counter;
	int matl_num;
	double rho, thickness, thickness_cubed, fdum, fdum2;
        double B_mass[MsoB], B2_mass[MsoB];
        double M_temp[neqlsq], M_el[neqlsq];
        double coord_el_trans[neqel];
        double det[num_int], area_el;
        double mass_el[neqel];

/*      initialize all variables  */

        memset(B_mass,0,MsoB*sof);
        memset(B2_mass,0,MsoB*sof);

	memcpy(shg.bend,shl.bend,soshb*sizeof(double));

        for( k = 0; k < numel; ++k )
        {
                matl_num = *(el_matl+k);
		thickness = matl[matl_num].thick;
		thickness_cubed = thickness*thickness*thickness/12.0;
                rho = matl[matl_num].rho;

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
			*(dof_el+ndof*j+2) = ndof*node+2;
                }

/* The call to plshg_mass is only for calculating the determinent */

		check = plshg_mass(det, k, shg, coord_el_trans, &area_el);
		if(!check) printf( "Problems with plshg_mass \n");

                /* printf("This is Area %10.6f for element %4d\n",area_el,k);*/

#if 0
                for( i1 = 0; i1 < num_int; ++i1 )
                {
                    for( i2 = 0; i2 < npel; ++i2 )
                    {
                    	printf("%10.6f ",*(shl+npel*(nsd+1)*i1 + npel*(nsd) + i2));
                    }
                    printf(" \n");
                }
                printf(" \n");
#endif

/* The loop over j below calculates the 4 points of the gaussian integration
   for several quantities */

                for( j = 0; j < num_int; ++j )
                {

/* Assembly of the B matrix for mass */

       		    check = plateB_mass((shg.bend+npel*(nsd+1)*j + npel*(nsd)),B_mass);
       		    if(!check) printf( "Problems with plateB_mass \n");

		    memcpy(B2_mass,B_mass,MsoB*sizeof(double));

                    check=matXT(M_temp, B_mass, B2_mass, neqel, neqel, nsd+1);
                    if(!check) printf( "Problems with matXT \n");

		    fdum = rho*thickness*(*(w+j))*(*(det+j));
		    fdum2 = rho*thickness_cubed*(*(w+j))*(*(det+j));

		    for( i2 = 0; i2 < npel; ++i2 )
		    {
			*(M_el+ndof*i2) +=
			    *(M_temp+ndof*i2)*fdum;
			*(M_el+neqel*1+ndof*i2+1) +=
			    *(M_temp+neqel*1+ndof*i2+1)*fdum2;
			*(M_el+neqel*2+ndof*i2+2) +=
			    *(M_temp+neqel*2+ndof*i2+2)*fdum2;

			*(M_el+neqel*3+ndof*i2) +=
			    *(M_temp+neqel*3+ndof*i2)*fdum;
			*(M_el+neqel*4+ndof*i2+1) +=
			    *(M_temp+neqel*4+ndof*i2+1)*fdum2;
			*(M_el+neqel*5+ndof*i2+2) +=
			    *(M_temp+neqel*5+ndof*i2+2)*fdum2;

			*(M_el+neqel*6+ndof*i2) +=
			    *(M_temp+neqel*6+ndof*i2)*fdum;
			*(M_el+neqel*7+ndof*i2+1) +=
			    *(M_temp+neqel*7+ndof*i2+1)*fdum2;
			*(M_el+neqel*8+ndof*i2+2) +=
			    *(M_temp+neqel*8+ndof*i2+2)*fdum2;

			*(M_el+neqel*9+ndof*i2) +=
			    *(M_temp+neqel*9+ndof*i2)*fdum;
			*(M_el+neqel*10+ndof*i2+1) +=
			    *(M_temp+neqel*10+ndof*i2+1)*fdum2;
			*(M_el+neqel*11+ndof*i2+2) +=
			    *(M_temp+neqel*11+ndof*i2+2)*fdum2;
		    }
                }

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
		    /*printf("This is Area %9.6f\n\n",fdum);*/

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
		/*printf("%5d  %16.8e\n",i, *(mass+i));*/
		if( *(id + i ) > -1 )
		{
		    *(mass + counter ) = *(mass + i );
		    ++counter;
		}
	    }
	}

        return 1;
}
