/*
    This program plots the mesh with the various
    forms of viewing including stress, strain, displacement
    materials, etc.  It works with a quad FEM code.

			San Le

 		Last Update 3/10/01

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
 */

#if WINDOWS
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../quad/qdconst.h"
#include "../quad/qdstruct.h"
#include "qdstrcgr.h"
#include "../../common_gr/control.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/* graphics globals */

extern choice2, choice3;

/* FEA globals */

extern double *coord, *coord0;
extern int *connecter;
extern int nmat, numnp, numel;
extern GLfloat MeshColor[boxnumber+5][3];
extern GLfloat wire_color[3], black[3], green[3], yellow[3];
extern GLfloat RenderColor[4];
extern ISTRESS *stress_color;
extern ISTRAIN *strain_color;
extern int *U_color, *el_matl_color;
extern int color_choice, input_flag, post_flag;
extern int input_color_flag;
extern int Solid_flag, Perspective_flag, Render_flag, AppliedDisp_flag,
        AppliedForce_flag, Material_flag, Node_flag, Element_flag, Axes_flag;
extern int Before_flag, After_flag, Both_flag, Amplify_flag;
extern int stress_flag, strain_flag, stress_strain, disp_flag;
extern int matl_choice, node_choice, ele_choice;

void qdmeshdraw(void)
{
        int i, i2, j, k, dof_el[neqel], sdof_el[npel*nsd], ii, check, counter, node;
	int l,m,n;
	int c0,c1,c2,c3;
	int matl_number, node_number;
	int After_gr_flag = 0, Before_gr_flag = 0;
        double coord_el[npel*3], coord0_el[npel*3], fpointx, fpointy, fpointz;
	GLfloat d1[3], d2[3], norm_temp[] = {0.0, 0.0, 1.0};

        if(post_flag + After_flag > 1) After_gr_flag = 1;
        else After_flag = 0;
        if(input_flag + Before_flag > 1) Before_gr_flag = 1;
        else Before_flag = 0;

	*(wire_color + 2) = 0.0;
	if(!Solid_flag) *(wire_color + 2) = 1.0;

        for( k = 0; k < numel; ++k )
        {
		for( j = 0; j < npel; ++j )
                {

/* Calculate element degrees of freedom */

                    node = *(connecter+npel*k+j);
                    *(sdof_el+nsd*j) = nsd*node;
                    *(sdof_el+nsd*j+1) = nsd*node+1;

		    *(dof_el+ndof*j) = ndof*node;
		    *(dof_el+ndof*j+1) = ndof*node+1;

/* Calculate local deformed coordinates */

		    if( post_flag )
		    {
                        *(coord_el+3*j)=*(coord+*(sdof_el+nsd*j));
                        *(coord_el+3*j+1)=*(coord+*(sdof_el+nsd*j+1));
                        *(coord_el+3*j+2)=0.0;
		    }

/* Calculate local undeformed coordinates */

		    if( input_flag )
		    {
                        *(coord0_el+3*j)=*(coord0+*(sdof_el+nsd*j));
                        *(coord0_el+3*j+1)=*(coord0+*(sdof_el+nsd*j+1));
                        *(coord0_el+3*j+2)=0.0;
		    }

    			/*printf( "%9.5f %9.5f %9.5f \n",*(coord_el+3*j),
				*(coord_el+3*j+1),*(coord_el+3*j+2));*/
                }

/* Calculate element material number */

		matl_number = *(el_matl_color + k);

        	switch (color_choice) {
               	    case 1:
			c0 = strain_color[k].pt[0].xx;
			c1 = strain_color[k].pt[1].xx;
			c2 = strain_color[k].pt[2].xx;
			c3 = strain_color[k].pt[3].xx;
               	    break;
               	    case 2:
			c0 = strain_color[k].pt[0].yy;
			c1 = strain_color[k].pt[1].yy;
			c2 = strain_color[k].pt[2].yy;
			c3 = strain_color[k].pt[3].yy;
               	    break;
               	    case 4:
			c0 = strain_color[k].pt[0].xy;
			c1 = strain_color[k].pt[1].xy;
			c2 = strain_color[k].pt[2].xy;
			c3 = strain_color[k].pt[3].xy;
               	    break;
               	    case 7:
			c0 = strain_color[k].pt[0].I;
			c1 = strain_color[k].pt[1].I;
			c2 = strain_color[k].pt[2].I;
			c3 = strain_color[k].pt[3].I;
               	    break;
               	    case 8:
			c0 = strain_color[k].pt[0].II;
			c1 = strain_color[k].pt[1].II;
			c2 = strain_color[k].pt[2].II;
			c3 = strain_color[k].pt[3].II;
               	    break;
               	    case 10:
			c0 = stress_color[k].pt[0].xx;
			c1 = stress_color[k].pt[1].xx;
			c2 = stress_color[k].pt[2].xx;
			c3 = stress_color[k].pt[3].xx;
               	    break;
               	    case 11:
			c0 = stress_color[k].pt[0].yy;
			c1 = stress_color[k].pt[1].yy;
			c2 = stress_color[k].pt[2].yy;
			c3 = stress_color[k].pt[3].yy;
               	    break;
               	    case 13:
			c0 = stress_color[k].pt[0].xy;
			c1 = stress_color[k].pt[1].xy;
			c2 = stress_color[k].pt[2].xy;
			c3 = stress_color[k].pt[3].xy;
               	    break;
               	    case 16:
			c0 = stress_color[k].pt[0].I;
			c1 = stress_color[k].pt[1].I;
			c2 = stress_color[k].pt[2].I;
			c3 = stress_color[k].pt[3].I;
               	    break;
               	    case 17:
			c0 = stress_color[k].pt[0].II;
			c1 = stress_color[k].pt[1].II;
			c2 = stress_color[k].pt[2].II;
			c3 = stress_color[k].pt[3].II;
               	    break;
               	    case 19:
			c0 = *(U_color + *(dof_el + ndof*0));
			c1 = *(U_color + *(dof_el + ndof*1));
			c2 = *(U_color + *(dof_el + ndof*2));
			c3 = *(U_color + *(dof_el + ndof*3));
               	    break;
               	    case 20:
			c0 = *(U_color + *(dof_el + ndof*0 + 1));
			c1 = *(U_color + *(dof_el + ndof*1 + 1));
			c2 = *(U_color + *(dof_el + ndof*2 + 1));
			c3 = *(U_color + *(dof_el + ndof*3 + 1));
               	    break;
               	    case 21:
			c0 = *(U_color + *(dof_el + ndof*0 + 1));
			c1 = *(U_color + *(dof_el + ndof*1 + 1));
			c2 = *(U_color + *(dof_el + ndof*2 + 1));
			c3 = *(U_color + *(dof_el + ndof*3 + 1));
               	    break;
               	    case 30:
			c0 = 0;
			c1 = 0;
			c2 = 0;
			c3 = 0;
			if( matl_choice == matl_number )
			{
				c0 = 7;
				c1 = 7;
				c2 = 7;
				c3 = 7;
			}
               	    break;
               	    case 31:
			c0 = 0;
			c1 = 0;
			c2 = 0;
			c3 = 0;
               	    break;
               	    case 32:
			c0 = 0;
			c1 = 0;
			c2 = 0;
			c3 = 0;
			if( ele_choice == k )
			{
				c0 = 7;
				c1 = 7;
				c2 = 7;
				c3 = 7;
			}
               	    break;
        	}

/* Draw the mesh after deformation */

    		if( After_gr_flag )
		{
    		   if( Solid_flag )
		   {
        	     	glBegin(GL_TRIANGLES);
                           	glNormal3fv(norm_temp);
			   	glColor3fv(MeshColor[c0]);
                	   	glVertex3dv((coord_el));
			   	glColor3fv(MeshColor[c3]);
                 	   	glVertex3dv((coord_el+9));
			   	glColor3fv(MeshColor[c1]);
                 	   	glVertex3dv((coord_el+3));
        	     	glEnd();
        	     	glBegin(GL_TRIANGLES);
                           	glNormal3fv(norm_temp);
			   	glColor3fv(MeshColor[c2]);
                 	   	glVertex3dv((coord_el+6));
			   	glColor3fv(MeshColor[c1]);
                 	   	glVertex3dv((coord_el+3));
			   	glColor3fv(MeshColor[c3]);
                 	   	glVertex3dv((coord_el+9));
        	     	glEnd();
		   }
   
/* Draw the wire frame around the mesh */
   
		   glColor3fv( black );
          	     glBegin(GL_LINE_LOOP);
                	   glVertex3dv((coord_el+9));
                 	   glVertex3dv((coord_el+6));
                 	   glVertex3dv((coord_el+3));
                 	   glVertex3dv((coord_el));
        	     glEnd();
		}

		if( input_color_flag )
		{
		     c0 = 8;
		     c1 = 8;
		     c2 = 8;
		     c3 = 8;
		}

/* Draw the mesh before deformation */

    		if( Before_gr_flag )
		{
    		   if( Solid_flag )
		   {
        	     	glBegin(GL_TRIANGLES);
                           	glNormal3fv(norm_temp);
			   	glColor3fv(MeshColor[c0]);
                	   	glVertex3dv((coord0_el));
			   	glColor3fv(MeshColor[c3]);
                 	   	glVertex3dv((coord0_el+9));
			   	glColor3fv(MeshColor[c1]);
                 	   	glVertex3dv((coord0_el+3));
        	     	glEnd();
        	     	glBegin(GL_TRIANGLES);
                           	glNormal3fv(norm_temp);
			   	glColor3fv(MeshColor[c2]);
                 	   	glVertex3dv((coord0_el+6));
			   	glColor3fv(MeshColor[c1]);
                 	   	glVertex3dv((coord0_el+3));
			   	glColor3fv(MeshColor[c3]);
                 	   	glVertex3dv((coord0_el+9));
        	     	glEnd();
		   }
   
/* Draw the wire frame around the mesh */
   
		   glColor3fv(wire_color);
          	     glBegin(GL_LINE_LOOP);
                	   glVertex3dv((coord0_el+9));
                 	   glVertex3dv((coord0_el+6));
                 	   glVertex3dv((coord0_el+3));
                 	   glVertex3dv((coord0_el));
        	     glEnd();
		}
	}
/* This draws the Node ID node */
	if (color_choice == 31)
	{
	    glPointSize(8);
	    node_number=node_choice;
    	    if( After_gr_flag )
	    {
	    	fpointx = *(coord+nsd*node_number);
	    	fpointy = *(coord+nsd*node_number+1);
	    	fpointz = 0.0;
              	glBegin(GL_POINTS);
        	    glColor3fv(yellow);
               	    glVertex3f(fpointx, fpointy, fpointz);
    	   	glEnd();
	    }
    	    if( Before_gr_flag )
	    {
	    	fpointx = *(coord0+nsd*node_number);
	    	fpointy = *(coord0+nsd*node_number+1);
	    	fpointz = 0.0;
              	glBegin(GL_POINTS);
        	    glColor3fv(yellow);
               	    glVertex3f(fpointx, fpointy, fpointz);
    	   	glEnd();
	    }
	}
	/*return 1;*/
}

