/*
    This utility function takes the product of a vector with the
    consistent mass matrix.  This is for a finite element program
    which does analysis on a beam.  It is for modal analysis.

		Updated 12/15/00

    SLFFEA source file
    Version:  1.3
    Copyright (C) 1999, 2000, 2001, 2002  San Le 

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

extern int numel, numnp, dof, sof;
extern int consistent_mass_flag, consistent_mass_store;

int bmnormcrossX(double *, double *, double *);

int matX(double *, double *, double *, int, int, int);

int matXrot(double *, double *, double *, int, int);

int rotXmat(double *, double *, double *, int, int);

int bmMassPassemble(double *axis_z, int *connect, double *coord, int *el_matl,
	int *el_type, double *mass, MATL *matl, double *P_global, double *U) 
{
        int i, i1, i2, i3, j, k, ij, dof_el[neqel];
        int check, counter, dum, node0, node1;
        int matl_num, el_num, type_num;
        double area, rho, Iydivarea, Izdivarea, Ipdivarea, fdum, fdum2;
        double L, Lx, Ly, Lz, Lsq, Lxysq, axis_x[nsd], axis_y[nsd];
        double M_temp[neqlsq], M_el[neqlsq], rotate[nsdsq], rotateT[nsdsq];
	double U_el[neqel];
        double jacob;
	double P_el[neqel];
        SHAPE sh;

	memset(P_global,0,dof*sof);

	if(consistent_mass_store)
	{

/* Assemble P matrix using stored element mass matrices */

	    for( k = 0; k < numel; ++k )
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
		type_num = *(el_type+k);
        	rho = matl[matl_num].rho;
        	area = matl[matl_num].area;
        	Iydivarea = matl[matl_num].Iy/area;
        	Izdivarea = matl[matl_num].Iz/area;
		Ipdivarea = (matl[matl_num].Iy + matl[matl_num].Iz)/area;

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
                if(!check)
		{

/* If magnitude of axis_y < SMALL(i.e. axis_z and axis_x are parallel), then make
local x global z, local z global x, and local y global y.  */

		   memset(rotate,0,nsdsq*sof);
		   memset(rotateT,0,nsdsq*sof);
                   *(rotate+2) = 1.0;
                   *(rotate+4) = 1.0;
                   *(rotate+6) = -1.0;
		   if(Lz < -ONE)
		   {
			memset(rotate,0,nsdsq*sof);
			memset(rotateT,0,nsdsq*sof);
                	*(rotate+2) = -1.0;
                	*(rotate+4) = 1.0;
                	*(rotate+6) = 1.0;
		   }
		   *(axis_z + nsd*k) = *(rotate+6);
		   *(axis_z + nsd*k+1) = 0.0;
		   *(axis_z + nsd*k+2) = 0.0;
		}
		else
		{

/* To find the true axis_z, take cross product of axis_x and axis_y */

		   check = bmnormcrossX(axis_x, axis_y, (axis_z+nsd*k));
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
		}

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

/* defining the components of a el element vector */

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

                memset(M_el,0,neqlsq*sof);

/* Assembly of the local mass matrix:

   For a truss, the only non-zero components are those for the axial DOF terms.  So
   leave as zero in [M_el] everything except *(M_el) and *(M_el+6)
   and *(M_el+72) and *(M_el+78) if the type_num = 1.

   Normally, we would form the product of [mass_el] = rho*[B_mass(transpose)][B_mass],
   but this cannot be done for beams because certain terms integrate to zero
   by analytical inspection, and so this matrix must be explicitly coded.  The zero
   terms pertain to integrations involving center of inertia, product of inertia, etc.
   The reason this can be done for the stiffness matrix is because the [B] matrix
   has a form which maintains the zeros of [K_el].

   Instead, I will use the fully integrated mass matrix as given by 
   "Theory of Matrix Structural Analysis" by J. S. Przemieniecki on page 295.
   
*/
/* row 0 */
		*(M_el) = 1.0/3.0;
		*(M_el+6) = 1.0/6.0;
		*(M_el+72) = *(M_el+6);
/* row 6 */
		*(M_el+78) = 1.0/3.0;

		if(type_num < 2)
		{
/* For truss element */

/* row 1 */
		    *(M_el+13) = 1.0/3.0;
		    *(M_el+19) = 1.0/6.0;
		    *(M_el+85) = *(M_el+19);
/* row 2 */
		    *(M_el+26) = 1.0/3.0;
		    *(M_el+32) = 1.0/6.0;
		    *(M_el+98) = *(M_el+32);
/* row 7 */
		    *(M_el+91) = 1.0/3.0;
/* row 8 */
		    *(M_el+104) = 1.0/3.0;
		}

		if(type_num > 1)
		{
/* row 1 */
			*(M_el+13) = 13.0/35.0 + 6.0*Izdivarea/(5.0*Lsq);

			*(M_el+17) = 11.0*L/210.0 + Izdivarea/(10.0*L);
			*(M_el+61) = *(M_el+17);

			*(M_el+19) = 9.0/70.0 - 6.0*Izdivarea/(5.0*Lsq);
			*(M_el+85) = *(M_el+19);

			*(M_el+23) =  -13.0*L/420.0 + Izdivarea/(10.0*L);
			*(M_el+133) = *(M_el+23);
/* row 2 */
			*(M_el+26) = 13.0/35.0 + 6.0*Iydivarea/(5.0*Lsq);

			*(M_el+28) = -11.0*L/210.0 - Iydivarea/(10.0*L);
			*(M_el+50) = *(M_el+28);

			*(M_el+32) = 9.0/70.0 - 6.0*Iydivarea/(5.0*Lsq);
			*(M_el+98) = *(M_el+32);

			*(M_el+34) = 13.0*L/420.0 - Iydivarea/(10.0*L);
			*(M_el+122) = *(M_el+34);
/* row 3 */
			*(M_el+39) = Ipdivarea/3.0;

			*(M_el+45) = Ipdivarea/6.0;
			*(M_el+111) = *(M_el+45);
/* row 4 */
			*(M_el+52) = Lsq/105.0 + 2.0*Iydivarea/15.0;

			*(M_el+56) = -13.0*L/420.0 + Iydivarea/(10.0*L);
			*(M_el+100) = *(M_el+56);

			*(M_el+58) = -Lsq/140.0 - Iydivarea/30.0;
			*(M_el+124) = *(M_el+58);
/* row 5 */
			*(M_el+65) = Lsq/105.0 + 2.0*Izdivarea/15.0;

			*(M_el+67) = 13.0*L/420.0 - Izdivarea/(10.0*L);
			*(M_el+89) = *(M_el+67);

			*(M_el+71) = -Lsq/140.0 - Izdivarea/30.0;
			*(M_el+137) = *(M_el+71);
/* row 7 */
			*(M_el+91) = 13.0/35.0 + 6.0*Izdivarea/(5.0*Lsq);

			*(M_el+95) = -11.0*L/210.0 - Izdivarea/(10.0*L);
			*(M_el+139) = *(M_el+95);
/* row 8 */
			*(M_el+104) = 13.0/35.0 + 6.0*Iydivarea/(5.0*Lsq);

			*(M_el+106) = 11.0*L/210.0 - Iydivarea/(10.0*L);
			*(M_el+128) = *(M_el+106);
/* row 9 */
			*(M_el+117) = Ipdivarea/3.0;
/* row 10 */
			*(M_el+130) =  Lsq/105.0 + 2.0*Iydivarea/15.0;
/* row 11 */
			*(M_el+143) =  Lsq/105.0 + 2.0*Izdivarea/15.0;
		}

		fdum = rho*area*L;
                for( i1 = 0; i1 < neqlsq; ++i1 )
                {
                        *(M_el + i1) = *(M_el + i1)*fdum;
                }

/* Put M_el back to global coordinates */

                check = matXrot(M_temp, M_el, rotate, neqel, neqel);
                if(!check) printf( "Problems with matXrot \n");

                check = rotXmat(M_el, rotateT, M_temp, neqel, neqel);
                if(!check) printf( "Problems with rotXmat \n");

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
