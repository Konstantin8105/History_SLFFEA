/*
    This program plots the rendered mesh for a plate.
  
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

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <stdlib.h>
#include "../plate/plconst.h"
#include "../plate/plstruct.h"
#include "plstrcgr.h"
#include "../../common_gr/control.h"

/* glut header files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/* graphics globals */

extern choice2, choice3;

/* FEA globals */

extern double *coord, *coord0, *zcoord;
extern NORM *norm, *norm0;
extern int *connecter;
extern int nmat, numnp, numel;
extern GLfloat MeshColor[boxnumber+5][3];
extern GLfloat RenderColor[4];
extern int color_flag, input_flag, post_flag;
extern int input_color_flag;
extern int Solid_flag, Perspective_flag, Render_flag, AppliedDisp_flag,
        AppliedForce_flag, Material_flag, Axes_flag;
extern int Before_flag, After_flag, Both_flag, Amplify_flag, Solid_flag,
 	stress_flag, strain_flag, stress_strain, disp_flag;

extern int matl_choice;

void plrender (void)
{
        int i, i2, j, k, sdof_el[npel*nsd], ii, check, counter, node;
	int l,m,n;
	int matl_number;
	int After_gr_flag = 0, Before_gr_flag = 0;
        double coord_el[npel*3], coord0_el[npel*3];
	GLfloat d1[3], d2[3], norm_temp[3];
        GLfloat frame_color;

	if(post_flag + After_flag > 1) After_gr_flag = 1;
	else After_flag = 0;
	if(input_flag + Before_flag > 1) Before_gr_flag = 1;
	else Before_flag = 0;
		
	frame_color = 0.0;
	if(!Solid_flag) frame_color = 1.0;

       	glMaterialfv(GL_FRONT, GL_DIFFUSE, RenderColor);

        for( k = 0; k < numel; ++k )
        {
                for( j = 0; j < npel; ++j )
                {

/* Calculate element degrees of freedom */

                        node = *(connecter+npel*k+j);
                        *(sdof_el+nsd*j) = nsd*node;
                        *(sdof_el+nsd*j+1) = nsd*node+1;

/* Calculate local deformed coordinates */

			if( post_flag )
			{
                        	*(coord_el+3*j)=*(coord+*(sdof_el+nsd*j));
                        	*(coord_el+3*j+1)=*(coord+*(sdof_el+nsd*j+1));
                        	*(coord_el+3*j+2)=*(zcoord+node);
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

/* Draw the mesh after deformation */

    		if( After_gr_flag )
		{
    		   if( Solid_flag )
		   {
			*(norm_temp) = norm[k].face[0].x;
			*(norm_temp+1) = norm[k].face[0].y;
			*(norm_temp+2) = norm[k].face[0].z;
        	        glBegin(GL_TRIANGLES);
                              glNormal3fv(norm_temp);
                	      glVertex3dv((coord_el));
                 	      glVertex3dv((coord_el+9));
                 	      glVertex3dv((coord_el+3));
        	        glEnd();
			*(norm_temp) = norm[k].face[1].x;
			*(norm_temp+1) = norm[k].face[1].y;
			*(norm_temp+2) = norm[k].face[1].z;
        	        glBegin(GL_TRIANGLES);
        	        glBegin(GL_TRIANGLES);
                              glNormal3fv(norm_temp);
                 	      glVertex3dv((coord_el+6));
                 	      glVertex3dv((coord_el+3));
                 	      glVertex3dv((coord_el+9));
        	        glEnd();
		   }
			      
/* Draw the wire frame around the mesh */
   
		   glColor3f( 0.0, 0.0, 0.0);
          	     glBegin(GL_LINE_LOOP);
                	   glVertex3dv((coord_el+9));
                 	   glVertex3dv((coord_el+6));
                 	   glVertex3dv((coord_el+3));
                 	   glVertex3dv((coord_el));
        	     glEnd();
		}

/* Draw the mesh before deformation */

    		if( Before_gr_flag )
		{
    		   if( Solid_flag )
		   {
			*(norm_temp) = norm0[k].face[0].x;
			*(norm_temp+1) = norm0[k].face[0].y;
			*(norm_temp+2) = norm0[k].face[0].z;
        	        glBegin(GL_TRIANGLES);
                              glNormal3fv(norm_temp);
                	      glVertex3dv((coord0_el));
                 	      glVertex3dv((coord0_el+9));
                 	      glVertex3dv((coord0_el+3));
        	        glEnd();
			*(norm_temp) = norm0[k].face[1].x;
			*(norm_temp+1) = norm0[k].face[1].y;
			*(norm_temp+2) = norm0[k].face[1].z;
        	        glBegin(GL_TRIANGLES);
                              glNormal3fv(norm_temp);
                 	      glVertex3dv((coord0_el+6));
                 	      glVertex3dv((coord0_el+3));
                 	      glVertex3dv((coord0_el+9));
        	        glEnd();
		   }
   
/* Draw the wire frame around the mesh */
   
		   glColor3f( 0.0, 0.0, frame_color);
          	     glBegin(GL_LINE_LOOP);
                	   glVertex3dv((coord0_el+9));
                 	   glVertex3dv((coord0_el+6));
                 	   glVertex3dv((coord0_el+3));
                 	   glVertex3dv((coord0_el));
        	     glEnd();
		}

/* This plots the structure */
	}
	/*return 1;*/
}

