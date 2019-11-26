/*
    This utility code uses the conjugate gradient method
    to solve the linear system [K][U] = [f] for a finite
    element program which does analysis on a beam.  The first function
    assembles the P matrix.  It is called by the second function
    which allocates the memory and goes through the steps of the algorithm.
    These go with the calculation of displacement.

		Updated 12/11/02

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
#include "bmconst.h"
#include "bmstruct.h"
#include "bmshape_struct.h"

#define SMALL      1.e-20

extern int analysis_flag, dof, numel, numnp, sof;
extern int lin_algebra_flag, numel_K, numel_P;
extern double x[num_int], x_node[num_int], w[num_int];
extern int  iteration_max, iteration_const, iteration;
extern double tolerance;

int bmnormcrossX(double *, double *, double *);

int bmshape(SHAPE *, double , double , double );

int matX(double *, double *, double *, int, int, int);

int matXT(double *, double *, double *, int, int, int);

int matXrot(double *, double *, double *, int, int);

int rotXmat(double *, double *, double *, int, int);

int dotX(double *, double *, double *, int);

int bmBoundary( double *, BOUND );


int bmConjPassemble(double *A, double *axis_z, int *connect, double *coord, int *el_matl,
	int *el_type, MATL *matl, double *P_global, double *U)
{
/* This function assembles the P_global matrix for the displacement calculation by
   taking the product [K_el]*[U_el].  Some of the [K_el] is stored in [A].

			Updated 7/8/00
*/
	int i, i1, i2, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, node0, node1;
	int matl_num, type_num;
	double area, Emod, EmodXarea, EmodXIy, EmodXIz, G, GXIp;
	double L, Lx, Ly, Lz, Lsq, Lxysq, axis_x[nsd], axis_y[nsd];
	double B[soB], DB[soB], jacob;
	double K_temp[neqlsq], K_el[neqlsq], K_local[neqlsq],
		rotate[nsdsq], rotateT[nsdsq];
	SHAPE sh;
	double U_el[neqel];
	double coord_el_trans[npel*nsd];
	double P_el[neqel];


	memset(P_global,0,dof*sof);

	for( k = 0; k < numel_K; ++k )
	{
		node0 = *(connect+k*npel);
		node1 = *(connect+k*npel+1);

		*(dof_el) = ndof*node0;
		*(dof_el+1) = ndof*node0+1;
		*(dof_el+2) = ndof*node0+2;
		*(dof_el+3) = ndof*node0+3;
		*(dof_el+4) = ndof*node0+4;
		*(dof_el+5) = ndof*node0+5;

		*(dof_el+6) = ndof*node1;
		*(dof_el+7) = ndof*node1+1;
		*(dof_el+8) = ndof*node1+2;
		*(dof_el+9) = ndof*node1+3;
		*(dof_el+10) = ndof*node1+4;
		*(dof_el+11) = ndof*node1+5;

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
		type_num = *(el_type+k);
		Emod = matl[matl_num].E;
		area = matl[matl_num].area;
		EmodXarea = matl[matl_num].E*matl[matl_num].area;
		EmodXIy = matl[matl_num].E*matl[matl_num].Iy;
		EmodXIz = matl[matl_num].E*matl[matl_num].Iz;
		G = matl[matl_num].E/(1.0 + matl[matl_num].nu)/2.0;
		GXIp = G*(matl[matl_num].Iy + matl[matl_num].Iz);

		node0 = *(connect+k*npel);
		node1 = *(connect+k*npel+1);
		Lx = *(coord+nsd*node1) - *(coord+nsd*node0);
		Ly = *(coord+nsd*node1+1) - *(coord+nsd*node0+1);
		Lz = *(coord+nsd*node1+2) - *(coord+nsd*node0+2);

		/*printf(" Lx, Ly, Lz %f %f %f\n ", Lx, Ly, Lz);*/

		Lsq = Lx*Lx+Ly*Ly+Lz*Lz;
		L = sqrt(Lsq);
		Lx /= L; Ly /= L; Lz /= L;
		*(axis_x) = Lx;
		*(axis_x+1) = Ly;
		*(axis_x+2) = Lz;

		jacob = L/2.0;

/* To find axis_y, take cross product of axis_z and axis_x */

		check = bmnormcrossX((axis_z+nsd*k), axis_x, axis_y);
		if(!check) printf( "Problems with bmnormcrossX \n");

/* Assembly of the 3X3 rotation matrix for the 12X12 global rotation
   matrix */
		memset(rotate,0,nsdsq*sof);
		*(rotate) = *(axis_x);
		*(rotate+1) = *(axis_x+1);
		*(rotate+2) = *(axis_x+2);
		*(rotate+3) = *(axis_y);
		*(rotate+4) = *(axis_y+1);
		*(rotate+5) = *(axis_y+2);
		*(rotate+6) = *(axis_z+nsd*k);
		*(rotate+7) = *(axis_z+nsd*k+1);
		*(rotate+8) = *(axis_z+nsd*k+2);


/* Assembly of the 3X3 transposed rotation matrix for the 12X12 global rotation
   matrix */

		memset(rotateT,0,nsdsq*sof);
		*(rotateT) = *(rotate);
		*(rotateT+1) = *(rotate+3);
		*(rotateT+2) = *(rotate+6);
		*(rotateT+3) = *(rotate+1);
		*(rotateT+4) = *(rotate+4);
		*(rotateT+5) = *(rotate+7);
		*(rotateT+6) = *(rotate+2);
		*(rotateT+7) = *(rotate+5);
		*(rotateT+8) = *(rotate+8);

/* defining the components of an el element vector */

		*(dof_el) = ndof*node0;
		*(dof_el+1) = ndof*node0+1;
		*(dof_el+2) = ndof*node0+2;
		*(dof_el+3) = ndof*node0+3;
		*(dof_el+4) = ndof*node0+4;
		*(dof_el+5) = ndof*node0+5;

		*(dof_el+6) = ndof*node1;
		*(dof_el+7) = ndof*node1+1;
		*(dof_el+8) = ndof*node1+2;
		*(dof_el+9) = ndof*node1+3;
		*(dof_el+10) = ndof*node1+4;
		*(dof_el+11) = ndof*node1+5;

		memset(U_el,0,neqel*sof);
		memset(K_el,0,neqlsq*sof);

/* The loop below calculates the 2 points of the gaussian integration */

		for( j = 0; j < num_int; ++j )
		{
		    memset(B,0,soB*sof);
		    memset(DB,0,soB*sof);
		    memset(K_local,0,neqlsq*sof);

		    check = bmshape(&sh, *(x+j), L, Lsq);
		    if(!check) printf( "Problems with bmshape \n");

/* Assembly of the local stiffness matrix.
   For a truss, the only non-zero components are those for
   the axial loads.  So leave as zero in [B] and [DB] everything except
   *(B) and *(B+6) and *(DB) and *(DB+6) if the type_num = 1.
*/

		    *(B) = sh.Nhat[0].dx1;
		    *(B+6) = sh.Nhat[1].dx1;
		    if(type_num > 1)
		    {
		    	*(B+13) = sh.N[0].dx2;
		    	*(B+17) = sh.N[1].dx2;
		    	*(B+19) = sh.N[2].dx2;
		    	*(B+23) = sh.N[3].dx2;
		    	*(B+26) = -sh.N[0].dx2;
		    	*(B+28) = sh.N[1].dx2;
		    	*(B+32) = -sh.N[2].dx2;
		    	*(B+34) = sh.N[3].dx2;
		    	*(B+39) = sh.Nhat[0].dx1;
		    	*(B+45) = sh.Nhat[1].dx1;
		    }

		    *(DB) = EmodXarea*sh.Nhat[0].dx1;
		    *(DB+6) = EmodXarea*sh.Nhat[1].dx1;
		    if(type_num > 1)
		    {
		    	*(DB+13) = EmodXIz*sh.N[0].dx2;
		    	*(DB+17) = EmodXIz*sh.N[1].dx2;
		    	*(DB+19) = EmodXIz*sh.N[2].dx2;
		    	*(DB+23) = EmodXIz*sh.N[3].dx2;
		    	*(DB+26) = -EmodXIy*sh.N[0].dx2;
		    	*(DB+28) = EmodXIy*sh.N[1].dx2;
		    	*(DB+32) = -EmodXIy*sh.N[2].dx2;
		    	*(DB+34) = EmodXIy*sh.N[3].dx2;
		    	*(DB+39) = GXIp*sh.Nhat[0].dx1;
		    	*(DB+45) = GXIp*sh.Nhat[1].dx1;
		    }

		    check = matXT(K_local, B, DB, neqel, neqel, sdim);
		    if(!check) printf( "Problems with matXT \n");

		    for( i1 = 0; i1 < neqlsq; ++i1 )
		    {
			*(K_el + i1) += *(K_local + i1)*jacob*(*(w+j));
		    }
		}

/* Put K back to global coordinates */

		check = matXrot(K_temp, K_el, rotate, neqel, neqel);
		if(!check) printf( "Problems with matXrot \n");

		check = rotXmat(K_el, rotateT, K_temp, neqel, neqel);
		if(!check) printf( "Problems with rotXmat \n");


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


int bmConjGrad(double *A, double *axis_z, BOUND bc, int *connect, double *coord,
	int *el_matl, int *el_type, double *force, double *K_diag, MATL *matl,
	double *U)
{
/* This function does memory allocation and uses the conjugate gradient
   method to solve the linear system arising from the calculation of
   displacements.  It also makes the call to bmConjPassemble to get the
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
	p=(mem_double+ptr_inc);	         ptr_inc += dof;
	P_global=(mem_double+ptr_inc);   ptr_inc += dof;
	rm1=(mem_double+ptr_inc);        ptr_inc += dof;
	r=(mem_double+ptr_inc);          ptr_inc += dof;
	z=(mem_double+ptr_inc);          ptr_inc += dof;
	zm1=(mem_double+ptr_inc);        ptr_inc += dof;

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
	check = bmBoundary (r, bc);
	if(!check) printf( " Problems with bmBoundary \n");

	check = bmBoundary (p, bc);
	if(!check) printf( " Problems with bmBoundary \n");

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
		check = bmConjPassemble( A, axis_z, connect, coord, el_matl, el_type,
			matl, P_global, p);
		if(!check) printf( " Problems with bmConjPassemble \n");
		check = bmBoundary (P_global, bc);
		if(!check) printf( " Problems with bmBoundary \n");
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
		check = bmBoundary (p, bc);
		if(!check) printf( " Problems with bmBoundary \n");

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
	check = bmConjPassemble( A, axis_z, connect, coord, el_matl, el_type, matl, P_global, U);
	if(!check) printf( " Problems with bmConjPassemble \n");

	for( j = 0; j < dof; ++j )
	{
		printf( "%4d %14.5f  %14.5f %14.5f  %14.5f  %14.5f %14.5f\n",j,alpha,beta,
			*(U+j),*(r+j),*(P_global+j),*(force+j));
	}
*/

	free(mem_double);

	return 1;
}


