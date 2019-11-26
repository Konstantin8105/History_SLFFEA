/*
    This utility function takes the product of a vector with the
    consistent mass matrix.  This is for a finite element program
    which does analysis on a truss.  It is for modal analysis.

		Updated 1/15/03

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
#include "tsconst.h"
#include "tsstruct.h"

extern int numel, numnp, dof, sof;
extern int consistent_mass_flag, consistent_mass_store;

int mat(double *, double *, double *, int, int, int);

int matXrot(double *, double *, double *, int, int);

int rotXmat(double *, double *, double *, int, int);

int tsMassPassemble(int *connect, double *coord, int *el_matl, double *mass,
	MATL *matl, double *P_global, double *U) 
{
        int i, i1, i2, i3, j, k, dof_el[neqel], sdof_el[npel*nsd];
	int check, counter, node0, node1;
	int matl_num;
	double area, rho, fdum, fdum2;
	double L, Lx, Ly, Lz, Lsq, Lxysq;
        double M_temp[neqlsq], M_el[neqlsq], rotate[nsdsq], rotateT[nsdsq];
	double U_el[neqel];
        double jacob;
	double P_el[neqel];

	memset(P_global,0,dof*sof);

	if(consistent_mass_store)
	{

/* Assemble P matrix using stored element mass matrices */

	    for( k = 0; k < numel; ++k )
	    {

		node0 = *(connect+k*npel);
		node1 = *(connect+k*npel+1);

		*(dof_el)=ndof*node0;
		*(dof_el+1)=ndof*node0+1;
		*(dof_el+2)=ndof*node0+2;
		*(dof_el+3)=ndof*node1;
		*(dof_el+4)=ndof*node1+1;
		*(dof_el+5)=ndof*node1+2;

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
        	area = matl[matl_num].area;

		node0 = *(connect+k*npel);
		node1 = *(connect+k*npel+1);
		Lx = *(coord+nsd*node1) - *(coord+nsd*node0);
		Ly = *(coord+nsd*node1+1) - *(coord+nsd*node0+1);
		Lz = *(coord+nsd*node1+2) - *(coord+nsd*node0+2);

                /*printf(" Lx, Ly, Lz %f %f %f\n ", Lx, Ly, Lz);*/

                Lsq = Lx*Lx+Ly*Ly+Lz*Lz;
                L = sqrt(Lsq);
		Lx /= L; Ly /= L; Lz /= L;
                jacob = L/2.0;

		Lxysq = Lx*Lx + Ly*Ly;
		Lxysq = sqrt(Lxysq);

/* Assembly of the 3X3 rotation matrix for the 6X6 global rotation
   matrix.  The mechanism below for calculating rotations is based on the method
   givin in the book, "A First Course in the Finite Element
   Method 2nd Ed." by Daryl L. Logan and my own.  See pages 236-239 in:

     Logan, Daryl L., A First Course in the Finite Element Method 2nd Ed., PWS-KENT,
        1992.
 */

		memset(rotate,0,nsdsq*sof);
		*(rotate) = Lx;
		*(rotate+1) = Ly;
		*(rotate+2) = Lz;
		*(rotate+3) = -Ly/Lxysq;
		*(rotate+4) = Lx/Lxysq;
		*(rotate+5) = 0.0;
		*(rotate+6) = -Lx*Lz/Lxysq;
		*(rotate+7) = -Ly*Lz/Lxysq;
		*(rotate+8) = Lxysq;

/* Assembly of the 3X3 transposed rotation matrix for the 6X6 global rotation
   matrix */

		memset(rotateT,0,nsdsq*sof);
		*(rotateT) = Lx;
		*(rotateT+1) = -Ly/Lxysq;
		*(rotateT+2) = -Lx*Lz/Lxysq;
		*(rotateT+3) = Ly;
		*(rotateT+4) = Lx/Lxysq;
		*(rotateT+5) = -Ly*Lz/Lxysq;
		*(rotateT+6) = Lz;
		*(rotateT+7) = 0.0;
		*(rotateT+8) = Lxysq;

/* If local x axis lies on global z axis, modify rotate and rotateT */

		if(Lz > ONE)
		{
			/*printf("%d %f\n", k, *(rotate+2)); */
			memset(rotate,0,nsdsq*sof);
			memset(rotateT,0,nsdsq*sof);
			*(rotate+2) = 1.0;
			*(rotate+4) = 1.0;
			*(rotate+6) = -1.0;

			*(rotateT+2) = -1.0;
			*(rotateT+4) = 1.0;
			*(rotateT+6) = 1.0;
		}

		if(Lz < -ONE)
		{
			/*printf("%d %f\n", k, *(rotate+2)); */
			memset(rotate,0,nsdsq*sof);
			memset(rotateT,0,nsdsq*sof);
			*(rotate+2) = -1.0;
			*(rotate+4) = 1.0;
			*(rotate+6) = 1.0;

			*(rotateT+2) = 1.0;
			*(rotateT+4) = 1.0;
			*(rotateT+6) = -1.0;
		}

/* defining the components of a el element vector */

		*(dof_el)=ndof*node0;
		*(dof_el+1)=ndof*node0+1;
		*(dof_el+2)=ndof*node0+2;
		*(dof_el+3)=ndof*node1;
		*(dof_el+4)=ndof*node1+1;
		*(dof_el+5)=ndof*node1+2;

                memset(M_el,0,neqlsq*sof);

		fdum = .5*rho*area*L;

/* The mass matrix below comes from equation 10.81 on page 279 of:

     Przemieniecki, J. S., Theory of Matrix Structural Analysis, Dover
        Publications Inc., New York, 1985.
*/

/* row 0 */
		*(M_el) = 1.0/3.0;
		*(M_el+3) = 1.0/6.0;
		*(M_el+18) = *(M_el+3);
/* row 1 */
		*(M_el+7) = 1.0/3.0;
		*(M_el+10) = 1.0/6.0;
		*(M_el+25) = *(M_el+10);
/* row 2 */
		*(M_el+14) = 1.0/3.0;
		*(M_el+17) = 1.0/6.0;
		*(M_el+32) = *(M_el+17);
/* row 3 */
		*(M_el+21) = 1.0/3.0;
/* row 4 */
		*(M_el+28) = 1.0/3.0;
/* row 5 */
		*(M_el+35) = 1.0/3.0;

		fdum = rho*area*L;
		for( i1 = 0; i1 < neqlsq; ++i1 )
		{
			*(M_el + i1) = *(M_el + i1)*fdum;
		}

/* Put M_el back to global coordinates */

		check = matXrot(M_temp, M_el, rotate, neqel, neqel);
		if(!check) printf( "error with matXrot \n");

		check = rotXmat(M_el, rotateT, M_temp, neqel, neqel);
		if(!check) printf( "error with rotXmat \n");

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
