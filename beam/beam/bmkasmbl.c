/*
    This utility function assembles the stiffness matrix for a finite 
    element program which does analysis on a beam.

		 Updated 10/27/00

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
#include "bmconst.h"
#include "bmstruct.h"
#include "bmshape_struct.h"

extern int analysis_flag, dof, numel, numnp, neqn, sof;
extern int gauss_stress_flag;
extern int lin_algebra_flag, numel_K, numel_P, numel_surf;
extern double x[num_int], x_node[num_int], w[num_int];

int globalConjKassemble(double *, int *, int , double *,
        double *, int , int , int );

int globalKassemble(double *, int *, double *, int *, int );

int matXrot(double *, double *, double *, int, int);

int rotXmat(double *, double *, double *, int, int);

int matX(double *, double *, double *, int, int, int);

int matXT(double *, double *, double *, int, int, int);

int bmshape(SHAPE *, double , double , double );

int bmnormcrossX(double *, double *, double *);

int bmKassemble(double *A, double *axis_z, BOUND bc, int *connect, double *coord,
	CURVATURE *curve, double *dist_load, int *el_matl, int *el_type,
	double *force, int *id, int *idiag, double *K_diag, int *lm, MATL *matl,
	MOMENT *moment, STRAIN *strain, STRESS *stress, double *U)
{
        int i, i1, j, k, ij, dof_el[neqel];
	int check, counter, dum, dist_flag, node, node0, node1;
        int matl_num, el_num, type_num;
        double area, Emod, EmodXarea, EmodXIy, EmodXIz, G, GXIp;
        double L, Lx, Ly, Lz, Lsq, Lxysq, axis_x[nsd], axis_y[nsd];
        double B[soB], DB[soB], jacob;
        double K_temp[neqlsq], K_el[neqlsq], K_local[neqlsq],
		rotate[nsdsq], rotateT[nsdsq];
        double force_el_dist[neqel], force_gl_dist[neqel],
		force_el_U[neqel], U_el[neqel], U_local[neqel];
	double stress_el[sdim], strain_el[sdim];
	SHAPE sh, sh_node;

        dum = 0;
        for( k = 0; k < numel; ++k )
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

/* The mechanism below for calculating rotations is based on a combination
   of the method givin in the book, "A First Course in the Finite Element
   Method 2nd Ed." by Daryl L. Logan and my own.  See pages 236-239 in:

     Logan, Daryl L., A First Course in the Finite Element Method 2nd Ed., PWS-KENT,
        1992.
*/

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
                   *(rotate+2) = 1.0;
                   *(rotate+4) = 1.0;
                   *(rotate+6) = -1.0;
		   if(Lz < -ONE)
		   {
			memset(rotate,0,nsdsq*sof);
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
                memset(force_el_dist,0,neqel*sof);
                memset(force_gl_dist,0,neqel*sof);
                memset(force_el_U,0,neqel*sof);

/* The loop below calculates the 2 points of the gaussian integration
   for several quantities */

                for( j = 0; j < num_int; ++j )
                {
                    memset(B,0,soB*sof);
                    memset(DB,0,soB*sof);
                    memset(K_local,0,neqlsq*sof);

                    check = bmshape(&sh, *(x+j), L, Lsq);
                    if(!check) printf( "Problems with bmshape \n");

/* Assembly of the local stiffness matrix:

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

/* assemble the distributed load */

                    if(analysis_flag == 1)
                    {
                      	if( k == bc.dist_load[dum] )
                      	{
		      	    el_num = bc.dist_load[dum];
                      	    *(force_el_dist) += 0.0;
                      	    *(force_el_dist+1) +=
				sh.N[0].dx0*(*(dist_load+2*el_num))*jacob*(*(w+j));
                      	    *(force_el_dist+2) +=
				sh.N[0].dx0*(*(dist_load+2*el_num+1))*jacob*(*(w+j));
                      	    *(force_el_dist+3) += 0.0;
                      	    *(force_el_dist+4) -=
				sh.N[1].dx0*(*(dist_load+2*el_num+1))*jacob*(*(w+j));
                      	    *(force_el_dist+5) +=
				sh.N[1].dx0*(*(dist_load+2*el_num))*jacob*(*(w+j));
                      	    *(force_el_dist+6) += 0.0;
                      	    *(force_el_dist+7) +=
				sh.N[2].dx0*(*(dist_load+2*el_num))*jacob*(*(w+j));
                      	    *(force_el_dist+8) +=
				sh.N[2].dx0*(*(dist_load+2*el_num+1))*jacob*(*(w+j));
                      	    *(force_el_dist+9) += 0.0;
                      	    *(force_el_dist+10) -=
				sh.N[3].dx0*(*(dist_load+2*el_num+1))*jacob*(*(w+j));
                      	    *(force_el_dist+11) +=
				sh.N[3].dx0*(*(dist_load+2*el_num))*jacob*(*(w+j));
                            /*printf ("bc. dist_load %d %d %d %f \n", i, el_num, bc.dist_load[dum],
                                *(dist_load+2*dum));*/
                      	}
                    }

		}
		dist_flag = 0;
                if( k == bc.dist_load[dum] )
		{
			++dum;
			dist_flag = 1;
		}

/* Put K back to global coordinates */

                check = matXrot(K_temp, K_el, rotate, neqel, neqel);
                if(!check) printf( "Problems with matXrot \n");

                check = rotXmat(K_el, rotateT, K_temp, neqel, neqel);
                if(!check) printf( "Problems with rotXmat \n");

                for( j = 0; j < neqel; ++j )
                {
                        *(U_el + j) = *(U + *(dof_el+j));
                }

                check = matX(force_el_U, K_el, U_el, neqel, 1, neqel);
                if(!check) printf( "Problems with matX \n");

                if(analysis_flag == 1)
                {

/* Compute the equivalant nodal forces based on prescribed displacements */

		  for( j = 0; j < neqel; ++j )
		  {
		  	*(force + *(dof_el+j)) -= *(force_el_U + j);
		  }

/* Compute the equivalant nodal forces based on distributed element loads */

		  if( dist_flag )
		  {
		     check = rotXmat(force_gl_dist, rotateT, force_el_dist, 1, neqel);
		     for( j = 0; j < neqel; ++j )
		     {
		     	*(force + *(dof_el+j)) += *(force_gl_dist + j);
		     }
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
			        *(force + *(dof_el+j)) += *(force_el_U + j);
			}

/* Calculate the element local U matrix */

                	check = rotXmat(U_local, rotate, U_el, 1, neqel);
                	if(!check) printf( "Problems with rotXmat \n");

/* The loop below calculates either the points of 2 point gaussian integration
   or at the nodes for the [B] and [stress] matrices.
*/

                	for( j = 0; j < num_int; ++j )
                	{
                    		memset(B,0,soB*sof);
				memset(stress_el,0,sdim*sof);
                                memset(strain_el,0,sdim*sof);

				if(gauss_stress_flag)
				{
/* Calculate sh at integration point */
                    		    check = bmshape(&sh, *(x+j), L, Lsq);
                    		    if(!check) printf( "Problems with bmshape \n");
				}
				else
				{
/* Calculate sh at nodal point */
                    		    check = bmshape(&sh_node, *(x_node+j), L, Lsq);
                    		    if(!check) printf( "Problems with bmshape \n");
				}

/* Similarly for what was discussed above:
   For a truss, the only non-zero components are those for
   the axial loads.  So leave as zero in [B] everything except
   *(B) and *(B+6) if the type_num = 1.
*/
				if(gauss_stress_flag)
				{
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
				}
				else
				{
                    		    *(B) = sh_node.Nhat[0].dx1;
                    		    *(B+6) = sh_node.Nhat[1].dx1;
				    if(type_num > 1)
				    {
                    			*(B+13) = sh_node.N[0].dx2;
                    			*(B+17) = sh_node.N[1].dx2;
                    			*(B+19) = sh_node.N[2].dx2;
                    			*(B+23) = sh_node.N[3].dx2;
                    			*(B+26) = -sh_node.N[0].dx2;
                    			*(B+28) = sh_node.N[1].dx2;
                    			*(B+32) = -sh_node.N[2].dx2;
                    			*(B+34) = sh_node.N[3].dx2;
                    			*(B+39) = sh_node.Nhat[0].dx1;
                    			*(B+45) = sh_node.Nhat[1].dx1;
				    }
				}

/* Calculation of the local strain matrix */

				check=matX(strain_el, B, U_local, sdim, 1, neqel );
                        	if(!check) printf( "Problems with matX\n");

/* Update of the global strain matrix */

                        	strain[k].pt[j].xx = *(strain_el);
                        	curve[k].pt[j].zz = *(strain_el+1);
                        	curve[k].pt[j].yy = *(strain_el+2);
                        	curve[k].pt[j].xx = *(strain_el+3);

/* Calculation of the local stress matrix */

                        	*(stress_el) = strain[k].pt[j].xx*Emod;
                        	*(stress_el+1) = curve[k].pt[j].zz*EmodXIz;
                        	*(stress_el+2) = curve[k].pt[j].yy*EmodXIy;
                        	*(stress_el+3) = curve[k].pt[j].xx*GXIp;

/* Update of the global stress matrix */

                        	stress[k].pt[j].xx += *(stress_el);
                        	moment[k].pt[j].zz += *(stress_el+1);
                        	moment[k].pt[j].yy += *(stress_el+2);
                        	moment[k].pt[j].xx += *(stress_el+3);
			}

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

	return 1;
}

