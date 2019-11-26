/*
    This utility function assembles the Mass and force matrix for a finite 
    element program which does analysis on a brick which can behave
    nonlinearly.  The equations are solved using dynamic relaxation.   

		Updated 8/22/06

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../brick/brconst.h"
#include "../brick/brstruct.h"

extern int numel, numnp, dof, sof;
extern double shg[sosh], shgh[sosh], shl[sosh], w[num_int], *Vol0;

int matXT(double *, double *, double *, int, int, int);

int brickB(double *,double *);

int brshg( double *, int, double *, double *, double *);

int brFMassemble(int *connect, double *coord, double *coordh, int *el_matl,
	double *force, double *mass, MATL *matl, double *U) 
	
{
        int i, i1, i2, i3, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, node;
	int matl_num;
	double Emod, G, K, Pois;
        double lamda, mu;
        double B[soB], Bh[soB], DB[soB];
        double K_temp[neqlsq], K_el[neqlsq];
	double force_el[neqel], U_el[neqel];
        double coord_el_trans[neqel], coordh_el_trans[neqel];
        double det[num_int], deth[num_int], volume_el, wXdet;
        double mass_el[neqel];

/*      initialize all variables  */
        memset(mass,0,dof*sof);
        memset(B,0,soB*sof);
        memset(Bh,0,soB*sof);

        for( k = 0; k < numel; ++k )
        {
                matl_num = *(el_matl+k);
                Emod = matl[matl_num].E;
                Pois = matl[matl_num].nu;

                K=Emod/(1.0-2*Pois)/3.0;
                G=Emod/(1.0+Pois)/2.0;

		volume_el = 0.0;

        	lamda = Emod*Pois/((1.0+Pois)*(1.0-2.0*Pois));
        	mu = Emod/(1.0+Pois)/2.0;

/* Zero out the Element stiffness and mass matrices */

        	memset(K_el,0,neqlsq*sof);
        	memset(mass_el,0,neqel*sof);

                for( j = 0; j < npel; ++j )
                {
			node = *(connect+npel*k+j);

			*(sdof_el+nsd*j) = nsd*node;
			*(sdof_el+nsd*j+1) = nsd*node+1;
			*(sdof_el+nsd*j+2) = nsd*node+2;

/* Create the coord vector and coordh_trans for one element */

			*(coord_el_trans+j)=*(coord+*(sdof_el+nsd*j));
			*(coord_el_trans+npel*1+j)=*(coord+*(sdof_el+nsd*j+1));
			*(coord_el_trans+npel*2+j)=*(coord+*(sdof_el+nsd*j+2));

/* Create the coordh and coordh_trans vector for one element */

                        *(coordh_el_trans+j) = *(coordh+nsd*node);
                        *(coordh_el_trans+npel*1+j) = *(coordh+nsd*node+1);
                        *(coordh_el_trans+npel*2+j) = *(coordh+nsd*node+2);

			*(dof_el+ndof*j) = ndof*node;
			*(dof_el+ndof*j+1) = ndof*node+1;
			*(dof_el+ndof*j+2) = ndof*node+2;
                }

/* Assembly of the shg matrix for each integration point at 1/2 time */

		brshg(deth, k, shl, shgh, coordh_el_trans);

/* Assembly of the shgh matrix for each integration point at full time */

		brshg(det, k, shl, shg, coord_el_trans);

/* The loop over j below calculates the 8 points of the gaussian integration
   for several quantities */

                for( j = 0; j < num_int; ++j )
                {

/* Assembly of the Bh matrix at 1/2 time */

                    check = brickB((shgh+npel*(nsd+1)*j),Bh);
                    if(!check) printf( "Problems with brickB \n");

/* Assembly of the B matrix at full time */

                    check = brickB((shg+npel*(nsd+1)*j),B);
                    if(!check) printf( "Problems with brickB \n");

                    for( i1 = 0; i1 < neqel; ++i1 )
                    {
                        *(DB+i1) = *(Bh+i1)*(lamda+2.0*mu)+*(Bh+neqel*1+i1)*lamda+
                                *(Bh+neqel*2+i1)*lamda;
                        *(DB+neqel*1+i1) = *(Bh+i1)*lamda+*(Bh+neqel*1+i1)*(lamda+2.0*mu)+
                                *(Bh+neqel*2+i1)*lamda;
                        *(DB+neqel*2+i1) = *(Bh+i1)*lamda+*(Bh+neqel*1+i1)*lamda+
                                *(Bh+neqel*2+i1)*(lamda+2.0*mu);
                        *(DB+neqel*3+i1) = *(Bh+neqel*3+i1)*mu;
                        *(DB+neqel*4+i1) = *(Bh+neqel*4+i1)*mu;
                        *(DB+neqel*5+i1) = *(Bh+neqel*5+i1)*mu;
                    }

		    wXdet = *(w+j)*(*(det+j));

/* Calculate the Volume from determinant of the Jacobian */

		    volume_el += wXdet;

                    check=matXT(K_temp, B, DB, neqel, neqel, sdim);
                    if(!check) printf( "Problems with matXT\n");
                    for( i2 = 0; i2 < neqlsq; ++i2 )
                    {
                          *(K_el+i2) += *(K_temp+i2)*wXdet;
                    }
                }

		for( j = 0; j < neqel; ++j )
		{
			*(U_el + j) = *(U + *(dof_el+j));
		}

		check = matX(force_el, K_el, U_el, neqel, 1, neqel);
		if(!check) printf( "Problems with matX \n");

/* Compute the equivalant nodal forces based on prescribed displacements */

		for( j = 0; j < neqel; ++j )
		{
			*(force + *(dof_el+j)) -= *(force_el + j);
		}

/* Creating the mass Matrix */
                for( i3 = 0; i3 < neqel; ++i3 )
                {
		   *(mass_el+i3) = 100.0*(*(K_el+neqel*i3+i3));
                }
                for( j = 0; j < npel; ++j )
                {
                    *(mass+*(dof_el+ndof*j)) += *(mass_el + ndof*j);
                    *(mass+*(dof_el+ndof*j+1)) += *(mass_el + ndof*j + 1);
                    *(mass+*(dof_el+ndof*j+2)) += *(mass_el + ndof*j + 2);
		}
        }
	/*for( i = 0; i < dof ; ++i )
	{
		printf( " force %4d %16.4e \n",i,*(force+i));
	}*/
        return 1;
}
