/*
    This utility code uses the conjugate gradient method
    to solve the linear system [K][U] = [f] for a finite
    element program which does analysis on a plate.  The first function
    assembles the P matrix.  It is called by the second function
    which allocates the memory and goes through the steps of the algorithm.
    These go with the calculation of displacement.

		Updated 12/11/02

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

#define SMALL      1.e-20

extern int analysis_flag, dof, numel, numnp, sof;
extern int lin_algebra_flag, numel_K, numel_P;
extern SH shg, shl;
extern double w[num_int+1], *Area0;
extern int  iteration_max, iteration_const, iteration;
extern double tolerance;

int matX(double *,double *, double *, int ,int ,int );

int matXT(double *, double *, double *, int, int, int);

int plateB4pt( double *, double *);

int plateB4pt_node( double *, double *);

int plshg( double *, int, SH, SH, double *, double *);

int dotX(double *, double *, double *, int);

int plBoundary( double *, BOUND );


int plConjPassemble(double *A, int *connect, double *coord, int *el_matl, MATL *matl,
	double *P_global, double *U)
{
/* This function assembles the P_global matrix for the displacement calculation by
   taking the product [K_el]*[U_el].  Some of the [K_el] is stored in [A].

			Updated 7/3/00
*/
        int i, i1, i2, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, node;
	int matl_num;
	double Emod, Pois, G1, G2, G3, thickness, shearK, const1, const2;
	double fdum, fdum2;
	double D11,D12,D21,D22;
	double lamda, mu;
	double B[soB], DB[soB];
	double K_temp[neqlsq], K_el[neqlsq];
	double U_el[neqel];
        double coord_el_trans[npel*nsd];
	double det[num_int+1];
        double P_el[neqel];


	memset(P_global,0,dof*sof);

        for( k = 0; k < numel_K; ++k )
        {

                for( j = 0; j < npel; ++j )
                {
			node = *(connect+npel*k+j);

                	*(dof_el+ndof*j) = ndof*node;
                	*(dof_el+ndof*j+1) = ndof*node+1;
                	*(dof_el+ndof*j+2) = ndof*node+2;
		}


/* Assembly of the global P matrix */

		for( j = 0; j < neqel; ++j )
		{
			*(U_el + j) = *(U + *(dof_el+j));
		}

                check = matX(P_el, (A+k*neqlsq), U_el, neqel, 1, neqel);
                if(!check) printf( "Problems with matX \n");

                for( j = 0; j < neqel; ++j )
                {
                	*(P_global+*(dof_el+j)) += *(P_el+j);
		}
	}

        for( k = numel_K; k < numel; ++k )
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
		}


/* Assembly of the shg matrix for each integration point */

		check = plshg(det, k, shl, shg, coord_el_trans, &fdum);
		if(!check) printf( "Problems with plshg \n");

		memset(U_el,0,neqel*sof);
                memset(K_el,0,neqlsq*sof);

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

                    check = matXT(K_temp, B, DB, neqel, neqel, sdim);
                    if(!check) printf( "Problems with matXT \n");

/* Note that I'm using the 4X4 determinant and weight for the 1X1 Gauss integration
   which is added 4 times.  This is a valid operation which can be proven. */

                    for( i2 = 0; i2 < neqlsq; ++i2 )
                    {
                          *(K_el+i2) += *(K_temp+i2)*(*(w+j))*(*(det+j));
                    }
                }

		for( j = 0; j < neqel; ++j )
                {
			*(U_el + j) = *(U + *(dof_el+j));
		}

/* Assembly of the global P matrix */

		for( j = 0; j < neqel; ++j )
		{
			*(U_el + j) = *(U + *(dof_el+j));
		}

		check = matX(P_el, K_el, U_el, neqel, 1, neqel);
                if(!check) printf( "Problems with matX \n");

                for( j = 0; j < neqel; ++j )
                {
                	*(P_global+*(dof_el+j)) += *(P_el+j);
		}
        }

        return 1;
}


int plConjGrad(double *A, BOUND bc, int *connect, double *coord, int *el_matl,
	double *force, double *K_diag, MATL *matl, double *U)
{
/* This function does memory allocation and uses the conjugate gradient
   method to solve the linear system arising from the calculation of
   displacements.  It also makes the call to plConjPassemble to get the
   product of [A]*[p].

			Updated 12/11/02

   It is taken from the algorithm 10.3.1 given in "Matrix Computations",
   by Golub, page 534.
*/
	int i, j, sofmf, ptr_inc;
	int check, counter;
	double *mem_double;
	double *p, *P_global, *r, *rm1, *z, *zm1;
	double alpha, alpha2, beta;
	double fdum, fdum2;

/* For the doubles */
	sofmf = 6*dof;
	mem_double=(double *)calloc(sofmf,sizeof(double));

	if(!mem_double )
	{
		printf( "failed to allocate memory for doubles\n ");
		exit(1);
	}

/* For the Conjugate Gradient Method doubles */

	                                        ptr_inc = 0;
	p=(mem_double+ptr_inc);                 ptr_inc += dof;
	P_global=(mem_double+ptr_inc);          ptr_inc += dof;
	rm1=(mem_double+ptr_inc);               ptr_inc += dof;
	r=(mem_double+ptr_inc);                 ptr_inc += dof;
	z=(mem_double+ptr_inc);                 ptr_inc += dof;
	zm1=(mem_double+ptr_inc);               ptr_inc += dof;

/* Using Conjugate gradient method to find displacements */

	memset(P_global,0,dof*sof);
	memset(p,0,dof*sof);
	memset(r,0,dof*sof);
	memset(rm1,0,dof*sof);
	memset(z,0,dof*sof);
	memset(zm1,0,dof*sof);

        for( j = 0; j < dof; ++j )
	{
		*(K_diag + j) += SMALL;
		*(r+j) = *(force+j);
		*(z + j) = *(r + j)/(*(K_diag + j));
		*(p+j) = *(z+j);
	}
	check = plBoundary (r, bc);
	if(!check) printf( " Problems with plBoundary \n");

	check = plBoundary (p, bc);
	if(!check) printf( " Problems with plBoundary \n");

	alpha = 0.0;
	alpha2 = 0.0;
	beta = 0.0;
	fdum2 = 1000.0;
	counter = 0;
	check = dotX(&fdum, r, z, dof);

	printf("\n iteration %3d iteration max %3d \n", iteration, iteration_max);
        /*for( iteration = 0; iteration < iteration_max; ++iteration )*/
	while(fdum2 > tolerance && counter < iteration_max )
        {

		printf( "\n %3d %16.8e\n",counter, fdum2);
		check = plConjPassemble( A, connect, coord, el_matl, matl, P_global, p);
		if(!check) printf( " Problems with plConjPassemble \n");
		check = plBoundary (P_global, bc);
		if(!check) printf( " Problems with plBoundary \n");
		check = dotX(&alpha2, p, P_global, dof);	
		alpha = fdum/(SMALL + alpha2);

        	for( j = 0; j < dof; ++j )
		{
            	    /*printf( "%4d %14.5e  %14.5e  %14.5e  %14.5e  %14.5e %14.5e\n",j,alpha,
			beta,*(U+j),*(r+j),*(P_global+j),*(p+j));*/
		    *(rm1+j) = *(r + j); 
		    *(zm1+j) = *(z + j); 
    		    *(U+j) += alpha*(*(p+j));
		    *(r+j) -=  alpha*(*(P_global+j));
		    *(z + j) = *(r + j)/(*(K_diag + j));
		}

		check = dotX(&fdum2, r, z, dof);
		beta = fdum2/(SMALL + fdum);
		fdum = fdum2;
		
        	for( j = 0; j < dof; ++j )
        	{
       		    /*printf("\n  %3d %12.7f  %14.5f ",j,*(U+j),*(P_global+j));*/
            	    /*printf( "%4d %14.5f  %14.5f  %14.5f  %14.5f %14.5f\n",j,alpha,
			*(U+j),*(r+j),*(P_global+j),*(force+j));
            	    printf( "%4d %14.8f  %14.8f  %14.8f  %14.8f %14.8f\n",j,
			*(U+j)*bet,*(r+j)*bet,*(P_global+j)*alp/(*(mass+j)),
			*(force+j)*alp/(*(mass+j)));*/
            	    *(p+j) = *(z+j)+beta*(*(p+j));

		}
		check = plBoundary (p, bc);
		if(!check) printf( " Problems with plBoundary \n");

		++counter;
	}

	if(counter > iteration_max - 1 )
	{
		printf( "\nMaximum iterations %4d reached.  Residual is: %16.8e\n",counter,fdum2);
		printf( "Problem may not have converged during Conj. Grad.\n");
	}
/*
The lines below are for testing the quality of the calculation:

1) r should be 0.0
2) P_global( = A*U ) - force should be 0.0
*/

/*
	check = plConjPassemble( A, connect, coord, el_matl, matl, P_global, U);
	if(!check) printf( " Problems with plConjPassemble \n");

	for( j = 0; j < dof; ++j )
	{
		printf( "%4d %14.5f  %14.5f %14.5f  %14.5f  %14.5f %14.5f\n",j,alpha,beta,
			*(U+j),*(r+j),*(P_global+j),*(force+j));
	}
*/

	free(mem_double);

	return 1;
}


