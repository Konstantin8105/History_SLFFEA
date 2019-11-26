/*
    This utility function assembles the K matrix for a finite 
    element program which does analysis on a quad.

		Updated 9/27/01

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
#include "qdconst.h"
#include "qdstruct.h"

extern int analysis_flag, dof, neqn, numel, numnp, plane_stress_flag, sof;
extern int gauss_stress_flag;
extern int lin_algebra_flag, numel_K, numel_P;
extern double shg[sosh], shg_node[sosh], shl[sosh], shl_node[sosh],
	shl_node2[sosh_node2], w[num_int], *Area0;

int globalConjKassemble(double *, int *, int , double *,
	double *, int , int , int );

int globalKassemble(double *, int *, double *, int *, int );

int matX( double *,double *,double *, int ,int ,int );

int matXT( double *, double *, double *, int, int, int);

int quadB( double *,double *);

int qdshg( double *, int, double *, double *, double *);

int qdstress_shg( double *, int, double *, double *, double * );

int qdKassemble(double *A, int *connect, double *coord, int *el_matl, double *force,
	int *id, int *idiag, double *K_diag, int *lm, MATL *matl,
	double *node_counter, STRAIN *strain, SDIM *strain_node, STRESS *stress,
	SDIM *stress_node, double *U)
{
	int i, i1, i2, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, counter, node;
	int matl_num;
	double Emod, Pois, G;
	double D11,D12,D21,D22;
	double lamda, mu;
	double B[soB], DB[soB];
	double K_temp[neqlsq], K_el[neqlsq];
	double force_el[neqel], U_el[neqel];
	double coord_el_trans[npel*nsd];
	double stress_el[sdim], strain_el[sdim], xxaddyy, xxsubyy, xysq;
	double det[num_int], wXdet;

	for( k = 0; k < numel; ++k )
	{

		matl_num = *(el_matl+k);
		Emod = matl[matl_num].E;
		Pois = matl[matl_num].nu;

		mu = Emod/(1.0+Pois)/2.0;

/* The lamda below is for plane strain */

		lamda = Emod*Pois/((1.0+Pois)*(1.0-2.0*Pois));

/* Recalculate lamda for plane stress */

		if(plane_stress_flag)
			lamda = Emod*Pois/(1.0-Pois*Pois);

		/*printf("lamda, mu, Emod, Pois  %f %f %f %f \n", lamda, mu, Emod, Pois);*/

		D11 = lamda+2.0*mu;
		D12 = lamda;
		D21 = lamda;
		D22 = lamda+2.0*mu;

		G = mu;

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

/* Count the number of times a particular node is part of an element */

			if(analysis_flag == 1)
				*(node_counter + node) += 1.0;
		}


/* Assembly of the shg matrix for each integration point */

		check=qdshg(det, k, shl, shg, coord_el_trans);
		if(!check) printf( "Problems with qdshg \n");

/* The loop over j below calculates the 4 points of the gaussian integration 
   for several quantities */

		memset(U_el,0,neqel*sof);
		memset(K_el,0,neqlsq*sof);
		memset(force_el,0,neqel*sof);

		for( j = 0; j < num_int; ++j )
		{

		    memset(B,0,soB*sof);
		    memset(DB,0,soB*sof);
		    memset(K_temp,0,neqlsq*sof);

/* Assembly of the B matrix */

		    check = quadB((shg+npel*(nsd+1)*j),B);
		    if(!check) printf( "Problems with quadB \n");

		    for( i1 = 0; i1 < neqel; ++i1 )
		    {
			*(DB+i1) = *(B+i1)*D11+
				*(B+neqel*1+i1)*D12;
			*(DB+neqel*1+i1) = *(B+i1)*D21+
				*(B+neqel*1+i1)*D22;
			*(DB+neqel*2+i1) = *(B+neqel*2+i1)*G;
		    }

		    wXdet = *(w+j)*(*(det+j));

		    check=matXT(K_temp, B, DB, neqel, neqel, sdim);
		    if(!check) printf( "Problems with matXT  \n");
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

		if(analysis_flag == 1)
		{

/* Compute the equivalant nodal forces based on prescribed displacements */

			for( j = 0; j < neqel; ++j )
			{
				*(force + *(dof_el+j)) -= *(force_el + j);
			}

/* Assembly of either the global skylined stiffness matrix or numel_K of the
   element stiffness matrices if the Conjugate Gradient method is used */

			if(lin_algebra_flag)
			{
			    check = globalKassemble(A, idiag, K_el, (lm + k*neqel),
				neqel);
			    if(!check) printf( "Problems with globalKassemble \n");
			}
			else
			{
			    check = globalConjKassemble(A, dof_el, k, K_diag, K_el,
				neqel, neqlsq, numel_K);
			    if(!check) printf( "Problems with globalConjKassemble \n");
			}
		}
		else
		{
/* Calculate the element reaction forces */

			for( j = 0; j < neqel; ++j )
			{
				*(force + *(dof_el+j)) += *(force_el + j);
			}

/* Calculate the element stresses */

/* Assembly of the qdstress_shg matrix for each nodal point */

/*
   qdstress_shg calculates shg at the nodes.  It is more efficient than qdshg
   because it removes all the zero multiplications from the calculation of shg. 
   You can use either function when calculating shg at the nodes. 

			    check=qdstress_shg(det, k, shl_node2, shg_node, coord_el_trans);
			    check=qdshg(det, k, shl_node, shg_node, coord_el_trans);
*/

			if(gauss_stress_flag)
			{
/* Calculate shg at integration point */
			    check=qdshg(det, k, shl, shg, coord_el_trans);
			    if(!check) printf( "Problems with qdshg \n");
			}
			else
			{
/* Calculate shg at nodal point */
			    check=qdstress_shg(det, k, shl_node2, shg_node, coord_el_trans);
			    if(!check) printf( "Problems with qdstress_shg \n");
			}

			for( j = 0; j < num_int; ++j )
			{

			   memset(B,0,soB*sof);
			   memset(stress_el,0,sdim*sof);
			   memset(strain_el,0,sdim*sof);

			   node = *(connect+npel*k+j);

/* Assembly of the B matrix */

			   if(gauss_stress_flag)
			   {
/* Calculate B matrix at integration point */
				check = quadB((shg+npel*(nsd+1)*j),B);
				if(!check) printf( "Problems with quadB \n");
			   }
			   else
			   {
/* Calculate B matrix at nodal point */
				check = quadB((shg_node+npel*(nsd+1)*j),B);
				if(!check) printf( "Problems with quadB \n");
			   }
		
/* Calculation of the local strain matrix */

			   check=matX(strain_el, B, U_el, sdim, 1, neqel );
			   if(!check) printf( "Problems with matX \n");

#if 0
			   for( i1 = 0; i1 < sdim; ++i1 )
			   {
				printf("%12.8f ",*(stress_el+i1));
				/*printf("%12.2f ",*(stress_el+i1));
				printf("%12.8f ",*(B+i1));*/
			   }
			   printf("\n");
#endif

/* Update of the global strain matrix */

			   strain[k].pt[j].xx = *(strain_el);
			   strain[k].pt[j].yy = *(strain_el+1);
			   strain[k].pt[j].xy = *(strain_el+2);

/* Calculate the principal straines */

			   xxaddyy = .5*(strain[k].pt[j].xx + strain[k].pt[j].yy);
			   xxsubyy = .5*(strain[k].pt[j].xx - strain[k].pt[j].yy);
			   xysq = strain[k].pt[j].xy*strain[k].pt[j].xy;

			   strain[k].pt[j].I = xxaddyy + sqrt( xxsubyy*xxsubyy
				+ xysq);
			   strain[k].pt[j].II = xxaddyy - sqrt( xxsubyy*xxsubyy
				+ xysq);
			   /*printf("%14.6e %14.6e %14.6e\n",xxaddyy,xxsubyy,xysq);*/

/* Add all the strains for a particular node from all the elements which share that node */

			   strain_node[node].xx += strain[k].pt[j].xx;
			   strain_node[node].yy += strain[k].pt[j].yy;
			   strain_node[node].xy += strain[k].pt[j].xy;
			   strain_node[node].I += strain[k].pt[j].I;
			   strain_node[node].II += strain[k].pt[j].II;

/* Calculation of the local stress matrix */

			   *(stress_el) = strain[k].pt[j].xx*D11+
				strain[k].pt[j].yy*D12;
			   *(stress_el+1) = strain[k].pt[j].xx*D21+
				strain[k].pt[j].yy*D22;
			   *(stress_el+2) = strain[k].pt[j].xy*G;

/* Update of the global stress matrix */

			   stress[k].pt[j].xx += *(stress_el);
			   stress[k].pt[j].yy += *(stress_el+1);
			   stress[k].pt[j].xy += *(stress_el+2);

/* Calculate the principal stresses */

			   xxaddyy = .5*(stress[k].pt[j].xx + stress[k].pt[j].yy);
			   xxsubyy = .5*(stress[k].pt[j].xx - stress[k].pt[j].yy);
			   xysq = stress[k].pt[j].xy*stress[k].pt[j].xy;

			   stress[k].pt[j].I = xxaddyy + sqrt( xxsubyy*xxsubyy
				+ xysq);
			   stress[k].pt[j].II = xxaddyy - sqrt( xxsubyy*xxsubyy
				+ xysq);

/* Add all the stresses for a particular node from all the elements which share that node */

			   stress_node[node].xx += stress[k].pt[j].xx;
			   stress_node[node].yy += stress[k].pt[j].yy;
			   stress_node[node].xy += stress[k].pt[j].xy;
			   stress_node[node].I += stress[k].pt[j].I;
			   stress_node[node].II += stress[k].pt[j].II;

/*
			   printf("%14.6e ",stress[k].pt[j].xx);
			   printf("%14.6e ",stress[k].pt[j].yy);
			   printf("%14.6e ",stress[k].pt[j].xy);
			   printf( "\n");
*/
			}
			/*printf( "\n");*/
		}
	}

	if(analysis_flag == 1)
	{

/* Contract the global force matrix using the id array only if linear
   algebra is used. */

	  if(lin_algebra_flag)
	  {
	     counter = 0;
	     for( i = 0; i < dof ; ++i )
	     {
		if( *(id + i ) > -1 )
		{
			*(force + counter ) = *(force + i );
			++counter;
		}
	     }
	  }
	}
	if(analysis_flag == 2)
	{

/* Average all the stresses and strains at the nodes */

	  for( i = 0; i < numnp ; ++i )
	  {
		   strain_node[i].xx /= *(node_counter + i);
		   strain_node[i].yy /= *(node_counter + i);
		   strain_node[i].xy /= *(node_counter + i);
		   strain_node[i].I /= *(node_counter + i);
		   strain_node[i].II /= *(node_counter + i);

		   stress_node[i].xx /= *(node_counter + i);
		   stress_node[i].yy /= *(node_counter + i);
		   stress_node[i].xy /= *(node_counter + i);
		   stress_node[i].I /= *(node_counter + i);
		   stress_node[i].II /= *(node_counter + i);
	  }
	}

	return 1;
}

