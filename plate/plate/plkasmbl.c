/*
    This utility function assembles the K matrix for a finite 
    element program which does analysis on a plate.

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
#include "plconst.h"
#include "plstruct.h"

extern int analysis_flag, dof, neqn, numel, numnp, sof;
extern int gauss_stress_flag;
extern int lin_algebra_flag, numel_K, numel_P;
extern SH shg, shg_node, shl, shl_node;
extern double shl_node2[sosh_node2], w[num_int+1], *Area0;

int globalConjKassemble(double *, int *, int , double *,
	double *, int , int , int );

int globalKassemble(double *, int *, double *, int *, int );

int matX( double *,double *,double *, int ,int ,int );

int matXT( double *, double *, double *, int, int, int);

int plateB4pt( double *, double *);

int plateB4pt_node( double *, double *);

int plateB1pt( double *, double *);

int plshg( double *, int, SH, SH, double *);

int plstress_shg( double *, int, double *, double *, double * );

int plKassemble(double *A, int *connect, double *coord, CURVATURE *curve, MDIM *curve_node,
	int *el_matl, double *force, int *id, int *idiag, double *K_diag, int *lm,
	MATL *matl, MOMENT *moment, MDIM *moment_node, double *node_counter,
	STRAIN *strain, SDIM *strain_node, STRESS *stress, SDIM *stress_node, double *U)
{
        int i, i1, i2, i3, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, counter, node;
	int matl_num;
	double Emod, Pois, G1, G2, G3, thickness, shearK, const1, const2;
	double lamda, mu;
	double D11,D12,D21,D22;
	double B[soB], DB[soB];
	double K_temp[neqlsq], K_el[neqlsq];
	double force_el[neqel], U_el[neqel];
        double coord_el_trans[npel*nsd];
	double stress_el[sdim], strain_el[sdim], xxaddyy, xxsubyy, xysq;
	double det[num_int+1], wXdet;

        for( k = 0; k < numel; ++k )
        {

                matl_num = *(el_matl+k);
                Emod = matl[matl_num].E;
		Pois = matl[matl_num].nu;
                thickness = matl[matl_num].thick;
		shearK = matl[matl_num].shear;

		mu = Emod/(1.0+Pois)/2.0;

/* The lamda below is for plane stress */

		lamda = Emod*Pois/(1.0-Pois*Pois);

                const1 = thickness*thickness*thickness/12.0;
                const2 = thickness*shearK;

        	/*printf("lamda, mu, Emod, Pois %f %f %f %f \n", lamda, mu, Emod, Pois);*/

		D11 = (lamda+2.0*mu)*const1;
		D12 = lamda*const1;
		D21 = lamda*const1;
		D22 = (lamda+2.0*mu)*const1;
		G1 = mu*const1;
		G2 = mu*const2;
		G3 = mu*const2;

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

/* Count the number of times a particular node is part of an element */

			if(analysis_flag == 1)
				*(node_counter + node) += 1.0;
		}


/* Assembly of the shg matrix for each integration point */

		check=plshg(det, k, shl, shg, coord_el_trans);
		if(!check) printf( "Problems with plshg \n");

		memset(U_el,0,neqel*sof);
                memset(K_el,0,neqlsq*sof);
		memset(force_el,0,neqel*sof);

/* Calculate the lower part of the B and DB matrices for 1X1 point gaussian integration */

                memset((B+neqel*(sdim-2)),0,neqel*(sdim-3)*sof);
                memset((DB+neqel*(sdim-2)),0,neqel*(sdim-3)*sof);

                check = plateB1pt(shg.shear, (B+(sdim-2)*neqel));
                if(!check) printf( "Problems with plateB1pt \n");
                for( i1 = 0; i1 < neqel; ++i1 )
                {
                        *(DB+neqel*3+i1)=*(B+neqel*3+i1)*G2;
                        *(DB+neqel*4+i1)=*(B+neqel*4+i1)*G3;
                }

/* The loop over j below calculates the 4 points of the gaussian integration 
   for several quantities */

                for( j = 0; j < num_int; ++j )
                {
		    memset(B,0,neqel*(sdim-2)*sof);
		    memset(DB,0,neqel*(sdim-2)*sof);
                    memset(K_temp,0,neqlsq*sof);

/* Assembly of the B matrix */

		    check = plateB4pt((shg.bend+npel*(nsd+1)*j), B);
                    if(!check) printf( "Problems with plateB4pt \n");

                    for( i1 = 0; i1 < neqel; ++i1 )
                    {
                        *(DB+i1)=*(B+i1)*D11 +
                                *(B+neqel*1+i1)*D12;
                        *(DB+neqel*1+i1)=*(B+i1)*D21 +
                                *(B+neqel*1+i1)*D22;
                        *(DB+neqel*2+i1)=*(B+neqel*2+i1)*G1;
                    }

		    wXdet = *(w+j)*(*(det+j));

                    check=matXT(K_temp, B, DB, neqel, neqel, sdim);
                    if(!check) printf( "Problems with matXT \n");

/* Note that I'm using the 4X4 determinant and weight for the 1X1 Gauss integration
   which is added 4 times.  This is a valid operation which can be proven. */

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

                	const2 = shearK;
			G2 = mu*const2;
			G3 = mu*const2;

/*
   plstress_shg calculates shg at the nodes.  It is more efficient than plshg
   because it removes all the zero multiplications from the calculation of shg.
   You can use either function when calculating shg at the nodes.

			   check=plstress_shg(det, k, shl_node2, shg_node.bend, coord_el_trans);
			   check=plshg(det, k, shl_node, shg_node, coord_el_trans);
*/

			if(gauss_stress_flag)
			{
/* Calculate shg at integration point */
			   check=plshg(det, k, shl, shg, coord_el_trans);
			   if(!check) printf( "Problems with plshg \n");
			}
			else
			{
/* Calculate shg at nodal point */
			   check=plstress_shg(det, k, shl_node2, shg_node.bend, coord_el_trans);
			   if(!check) printf( "Problems with plshg \n");
			}

			memset((B+neqel*(sdim-2)),0,neqel*(sdim-3)*sof);

/* Assembly of the B matrix */

			check = plateB1pt(shg.shear, (B+(sdim-2)*neqel));
			if(!check) printf( "Problems with plateB1pt \n");

			for( j = 0; j < num_int; ++j )
			{
			   memset(B,0,neqel*(sdim-2)*sof);
			   memset(stress_el,0,sdim*sof);
			   memset(strain_el,0,sdim*sof);

			   node = *(connect+npel*k+j);

/* Assembly of the B matrix */

			   if(gauss_stress_flag)
			   {
				check = plateB4pt((shg.bend+npel*(nsd+1)*j), B);
				if(!check) printf( "Problems with plateB4pt \n");
			   }
			   else
			   {
/* Calculate the shear terms in B at nodes using shg_node */
				check = plateB4pt_node((shg_node.bend+npel*(nsd+1)*j), B);	
				if(!check) printf( "Problems with plateB4pt_node \n");
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

			   curve[k].pt[j].xx = *(strain_el);
			   curve[k].pt[j].yy = *(strain_el+1);
			   curve[k].pt[j].xy = *(strain_el+2);
			   strain[k].pt[j].zx = *(strain_el+3);
			   strain[k].pt[j].yz = *(strain_el+4);

/* Calculate the principal straines */

			   xxaddyy = .5*(curve[k].pt[j].xx + curve[k].pt[j].yy);
			   xxsubyy = .5*(curve[k].pt[j].xx - curve[k].pt[j].yy);
			   xysq = curve[k].pt[j].xy*curve[k].pt[j].xy;

			   curve[k].pt[j].I = xxaddyy + sqrt( xxsubyy*xxsubyy
				+ xysq);
			   curve[k].pt[j].II = xxaddyy - sqrt( xxsubyy*xxsubyy
				+ xysq);
			   /*printf("%14.6e %14.6e %14.6e\n",xxaddyy,xxsubyy,xysq);*/

/* Add all the curvatures and strains for a particular node from all the elements which
   share that node */

			   curve_node[node].xx += curve[k].pt[j].xx;
			   curve_node[node].yy += curve[k].pt[j].yy;
			   curve_node[node].xy += curve[k].pt[j].xy;
			   strain_node[node].zx += strain[k].pt[j].zx;
			   strain_node[node].yz += strain[k].pt[j].yz;
			   curve_node[node].I += curve[k].pt[j].I;
			   curve_node[node].II += curve[k].pt[j].II;

/* Calculation of the local stress matrix */

			   *(stress_el)=curve[k].pt[j].xx*D11 +
				curve[k].pt[j].yy*D12;
			   *(stress_el+1)=curve[k].pt[j].xx*D21 +
				curve[k].pt[j].yy*D22;
			   *(stress_el+2)=curve[k].pt[j].xy*G1;
			   *(stress_el+3)=strain[k].pt[j].zx*G2;
			   *(stress_el+4)=strain[k].pt[j].yz*G3;

/* Update of the global stress matrix */

			   moment[k].pt[j].xx += *(stress_el);
			   moment[k].pt[j].yy += *(stress_el+1);
			   moment[k].pt[j].xy += *(stress_el+2);
			   stress[k].pt[j].zx += *(stress_el+3);
			   stress[k].pt[j].yz += *(stress_el+4);

/* Calculate the principal stresses */

			   xxaddyy = .5*(moment[k].pt[j].xx + moment[k].pt[j].yy);
			   xxsubyy = .5*(moment[k].pt[j].xx - moment[k].pt[j].yy);
			   xysq = moment[k].pt[j].xy*moment[k].pt[j].xy;

			   moment[k].pt[j].I = xxaddyy + sqrt( xxsubyy*xxsubyy
				+ xysq);
			   moment[k].pt[j].II = xxaddyy - sqrt( xxsubyy*xxsubyy
				+ xysq);

/* Add all the moments and stresses for a particular node from all the elements which
   share that node */

			   moment_node[node].xx += moment[k].pt[j].xx;
			   moment_node[node].yy += moment[k].pt[j].yy;
			   moment_node[node].xy += moment[k].pt[j].xy;
			   stress_node[node].zx += stress[k].pt[j].zx;
			   stress_node[node].yz += stress[k].pt[j].yz;
			   moment_node[node].I += moment[k].pt[j].I;
			   moment_node[node].II += moment[k].pt[j].II;

/*
			   printf("%14.6e ",moment[k].pt[j].xx);
			   printf("%14.6e ",moment[k].pt[j].yy);
			   printf("%14.6e ",moment[k].pt[j].xy);
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

/* Average all the moments, stresses, curvatures, and strains at the nodes */

          for( i = 0; i < numnp ; ++i )
	  {
		   curve_node[i].xx /= *(node_counter + i);
		   curve_node[i].yy /= *(node_counter + i);
		   curve_node[i].xy /= *(node_counter + i);
		   strain_node[i].zx /= *(node_counter + i);
		   strain_node[i].yz /= *(node_counter + i);
		   curve_node[i].I /= *(node_counter + i);
		   curve_node[i].II /= *(node_counter + i);

		   moment_node[i].xx /= *(node_counter + i);
		   moment_node[i].yy /= *(node_counter + i);
		   moment_node[i].xy /= *(node_counter + i);
		   stress_node[i].zx /= *(node_counter + i);
		   stress_node[i].yz /= *(node_counter + i);
		   moment_node[i].I /= *(node_counter + i);
		   moment_node[i].II /= *(node_counter + i);
          }
        }

	return 1;
}

