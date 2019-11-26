/*
    This utility function assembles the lumped sum Mass matrix for a
    finite element program which does analysis on a truss.  It is for
    modal analysis.

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

extern int dof, numel, numnp, sof;
extern int consistent_mass_flag, consistent_mass_store, lumped_mass_flag;

int matX(double *, double *, double *, int, int, int);

int matXT(double *, double *, double *, int, int, int);

int matXrot(double *, double *, double *, int, int);

int rotXmat(double *, double *, double *, int, int);

int tsMassemble( int *connect, double *coord, int *el_matl, int *id, double *mass,
	MATL *matl)
{
        int i, i1, i2, i3, j, k, ij, dof_el[neqel];
	int check, counter, dum, node0, node1;
        int matl_num, el_num, type_num;
        double area, rho, fdum, fdum2;
        double L, Lx, Ly, Lz, Lsq, Lxysq;
        double M_temp[neqlsq], M_el[neqlsq], mass_el[neqel],
		mass_local[neqlsq], rotate[nsdsq], rotateT[nsdsq];
        double jacob;

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
                memset(mass_el,0,neqel*sof);

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
		    /*printf("This is Volume2 %9.6f\n\n",fdum);*/

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

