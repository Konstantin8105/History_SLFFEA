/*
    This program plots the mesh with the various
    forms of viewing including stress, strain, displacement
    materials, etc.  It works with a shell FEM code.

			San Le

 		Last Update 9/25/06

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
#include "../shell/shconst.h"
#include "../shell/shstruct.h"
#include "shstrcgr.h"
#include "../../common_gr/control.h"

/* glut header files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/* graphics globals */

extern choice2, choice3;

/* FEA globals */

extern double *coord, *coord0;
extern NORM *norm, *norm0;
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
extern int matl_choice, node_choice, ele_choice;

void shmeshdraw(void)
{
        int i, i2, j, k, dof_el[neqel], sdof_el[npel*nsd], ii, check, counter, node;
	int l,m,n;
	int c0,c1,c2,c3,c4,c5,c6,c7;
	int matl_number, node_number;
	int After_gr_flag = 0, Before_gr_flag = 0;
        double coord_el[npel*3], coord0_el[npel*3], fpointx, fpointy, fpointz;
	GLfloat d1[3], d2[3], norm_temp[3];

        if(post_flag + After_flag > 1) After_gr_flag = 1;
        else After_flag = 0;
        if(input_flag + Before_flag > 1) Before_gr_flag = 1;
        else Before_flag = 0;

	*(wire_color + 2) = 0.0;
	if(!Solid_flag) *(wire_color + 2) = 1.0;

        for( k = 0; k < numel; ++k )
        {
		for( j = 0; j < npell; ++j )
                {

/* Calculate element degrees of freedom */

                    node = *(connecter+npell*k+j);

                    *(sdof_el+nsd*j) = nsd*node;
                    *(sdof_el+nsd*j+1) = nsd*node+1;
                    *(sdof_el+nsd*j+2) = nsd*node+2;

		    *(sdof_el+nsd*npell+nsd*j) = nsd*(node+numnp);
		    *(sdof_el+nsd*npell+nsd*j+1) = nsd*(node+numnp)+1;
		    *(sdof_el+nsd*npell+nsd*j+2) = nsd*(node+numnp)+2;

		    *(dof_el+ndof*j) = ndof*node;
		    *(dof_el+ndof*j+1) = ndof*node+1;
		    *(dof_el+ndof*j+2) = ndof*node+2;
		    *(dof_el+ndof*j+3) = ndof*node+3;
		    *(dof_el+ndof*j+4) = ndof*node+4;

/* Calculate local deformed coordinates */

		    if( post_flag )
		    {
                        *(coord_el+3*j)=*(coord+*(sdof_el+nsd*j));
                        *(coord_el+3*j+1)=*(coord+*(sdof_el+nsd*j+1));
                        *(coord_el+3*j+2)=*(coord+*(sdof_el+nsd*j+2));

                        *(coord_el+3*npell+3*j)=
				*(coord+*(sdof_el+nsd*npell+nsd*j));
                        *(coord_el+3*npell+3*j+1)=
				*(coord+*(sdof_el+nsd*npell+nsd*j+1));
                        *(coord_el+3*npell+3*j+2)=
				*(coord+*(sdof_el+nsd*npell+nsd*j+2));
		    }

/* Calculate local undeformed coordinates */

		    if( input_flag )
		    {
                        *(coord0_el+3*j)=*(coord0+*(sdof_el+nsd*j));
                        *(coord0_el+3*j+1)=*(coord0+*(sdof_el+nsd*j+1));
                        *(coord0_el+3*j+2)=*(coord0+*(sdof_el+nsd*j+2));

                        *(coord0_el+3*npell+3*j)=
				*(coord0+*(sdof_el+nsd*npell+nsd*j));
                        *(coord0_el+3*npell+3*j+1)=
				*(coord0+*(sdof_el+nsd*npell+nsd*j+1));
                        *(coord0_el+3*npell+3*j+2)=
				*(coord0+*(sdof_el+nsd*npell+nsd*j+2));
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
			c4 = strain_color[k].pt[4].xx;
			c5 = strain_color[k].pt[5].xx;
			c6 = strain_color[k].pt[6].xx;
			c7 = strain_color[k].pt[7].xx;
               	    break;
               	    case 2:
			c0 = strain_color[k].pt[0].yy;
			c1 = strain_color[k].pt[1].yy;
			c2 = strain_color[k].pt[2].yy;
			c3 = strain_color[k].pt[3].yy;
			c4 = strain_color[k].pt[4].yy;
			c5 = strain_color[k].pt[5].yy;
			c6 = strain_color[k].pt[6].yy;
			c7 = strain_color[k].pt[7].yy;
               	    break;
               	    case 4:
			c0 = strain_color[k].pt[0].xy;
			c1 = strain_color[k].pt[1].xy;
			c2 = strain_color[k].pt[2].xy;
			c3 = strain_color[k].pt[3].xy;
			c4 = strain_color[k].pt[4].xy;
			c5 = strain_color[k].pt[5].xy;
			c6 = strain_color[k].pt[6].xy;
			c7 = strain_color[k].pt[7].xy;
               	    break;
               	    case 5:
			c0 = strain_color[k].pt[0].zx;
			c1 = strain_color[k].pt[1].zx;
			c2 = strain_color[k].pt[2].zx;
			c3 = strain_color[k].pt[3].zx;
			c4 = strain_color[k].pt[4].zx;
			c5 = strain_color[k].pt[5].zx;
			c6 = strain_color[k].pt[6].zx;
			c7 = strain_color[k].pt[7].zx;
               	    break;
               	    case 6:
			c0 = strain_color[k].pt[0].yz;
			c1 = strain_color[k].pt[1].yz;
			c2 = strain_color[k].pt[2].yz;
			c3 = strain_color[k].pt[3].yz;
			c4 = strain_color[k].pt[4].yz;
			c5 = strain_color[k].pt[5].yz;
			c6 = strain_color[k].pt[6].yz;
			c7 = strain_color[k].pt[7].yz;
               	    break;
               	    case 7:
			c0 = strain_color[k].pt[0].I;
			c1 = strain_color[k].pt[1].I;
			c2 = strain_color[k].pt[2].I;
			c3 = strain_color[k].pt[3].I;
			c4 = strain_color[k].pt[4].I;
			c5 = strain_color[k].pt[5].I;
			c6 = strain_color[k].pt[6].I;
			c7 = strain_color[k].pt[7].I;
               	    break;
               	    case 8:
			c0 = strain_color[k].pt[0].II;
			c1 = strain_color[k].pt[1].II;
			c2 = strain_color[k].pt[2].II;
			c3 = strain_color[k].pt[3].II;
			c4 = strain_color[k].pt[4].II;
			c5 = strain_color[k].pt[5].II;
			c6 = strain_color[k].pt[6].II;
			c7 = strain_color[k].pt[7].II;
               	    break;
               	    case 9:
			c0 = strain_color[k].pt[0].III;
			c1 = strain_color[k].pt[1].III;
			c2 = strain_color[k].pt[2].III;
			c3 = strain_color[k].pt[3].III;
			c4 = strain_color[k].pt[4].III;
			c5 = strain_color[k].pt[5].III;
			c6 = strain_color[k].pt[6].III;
			c7 = strain_color[k].pt[7].III;
               	    break;
               	    case 10:
			c0 = stress_color[k].pt[0].xx;
			c1 = stress_color[k].pt[1].xx;
			c2 = stress_color[k].pt[2].xx;
			c3 = stress_color[k].pt[3].xx;
			c4 = stress_color[k].pt[4].xx;
			c5 = stress_color[k].pt[5].xx;
			c6 = stress_color[k].pt[6].xx;
			c7 = stress_color[k].pt[7].xx;
               	    break;
               	    case 11:
			c0 = stress_color[k].pt[0].yy;
			c1 = stress_color[k].pt[1].yy;
			c2 = stress_color[k].pt[2].yy;
			c3 = stress_color[k].pt[3].yy;
			c4 = stress_color[k].pt[4].yy;
			c5 = stress_color[k].pt[5].yy;
			c6 = stress_color[k].pt[6].yy;
			c7 = stress_color[k].pt[7].yy;
               	    break;
               	    case 13:
			c0 = stress_color[k].pt[0].xy;
			c1 = stress_color[k].pt[1].xy;
			c2 = stress_color[k].pt[2].xy;
			c3 = stress_color[k].pt[3].xy;
			c4 = stress_color[k].pt[4].xy;
			c5 = stress_color[k].pt[5].xy;
			c6 = stress_color[k].pt[6].xy;
			c7 = stress_color[k].pt[7].xy;
               	    break;
               	    case 14:
			c0 = stress_color[k].pt[0].zx;
			c1 = stress_color[k].pt[1].zx;
			c2 = stress_color[k].pt[2].zx;
			c3 = stress_color[k].pt[3].zx;
			c4 = stress_color[k].pt[4].zx;
			c5 = stress_color[k].pt[5].zx;
			c6 = stress_color[k].pt[6].zx;
			c7 = stress_color[k].pt[7].zx;
               	    break;
               	    case 15:
			c0 = stress_color[k].pt[0].yz;
			c1 = stress_color[k].pt[1].yz;
			c2 = stress_color[k].pt[2].yz;
			c3 = stress_color[k].pt[3].yz;
			c4 = stress_color[k].pt[4].yz;
			c5 = stress_color[k].pt[5].yz;
			c6 = stress_color[k].pt[6].yz;
			c7 = stress_color[k].pt[7].yz;
               	    break;
               	    case 16:
			c0 = stress_color[k].pt[0].I;
			c1 = stress_color[k].pt[1].I;
			c2 = stress_color[k].pt[2].I;
			c3 = stress_color[k].pt[3].I;
			c4 = stress_color[k].pt[4].I;
			c5 = stress_color[k].pt[5].I;
			c6 = stress_color[k].pt[6].I;
			c7 = stress_color[k].pt[7].I;
               	    break;
               	    case 17:
			c0 = stress_color[k].pt[0].II;
			c1 = stress_color[k].pt[1].II;
			c2 = stress_color[k].pt[2].II;
			c3 = stress_color[k].pt[3].II;
			c4 = stress_color[k].pt[4].II;
			c5 = stress_color[k].pt[5].II;
			c6 = stress_color[k].pt[6].II;
			c7 = stress_color[k].pt[7].II;
               	    break;
               	    case 18:
			c0 = stress_color[k].pt[0].III;
			c1 = stress_color[k].pt[1].III;
			c2 = stress_color[k].pt[2].III;
			c3 = stress_color[k].pt[3].III;
			c4 = stress_color[k].pt[4].III;
			c5 = stress_color[k].pt[5].III;
			c6 = stress_color[k].pt[6].III;
			c7 = stress_color[k].pt[7].III;
               	    break;
               	    case 19:
			c0 = *(U_color + *(dof_el + ndof*0));
			c1 = *(U_color + *(dof_el + ndof*1));
			c2 = *(U_color + *(dof_el + ndof*2));
			c3 = *(U_color + *(dof_el + ndof*3));
			c4 = *(U_color + *(dof_el + ndof*0));
			c5 = *(U_color + *(dof_el + ndof*1));
			c6 = *(U_color + *(dof_el + ndof*2));
			c7 = *(U_color + *(dof_el + ndof*3));
               	    break;
               	    case 20:
			c0 = *(U_color + *(dof_el + ndof*0 + 1));
			c1 = *(U_color + *(dof_el + ndof*1 + 1));
			c2 = *(U_color + *(dof_el + ndof*2 + 1));
			c3 = *(U_color + *(dof_el + ndof*3 + 1));
			c4 = *(U_color + *(dof_el + ndof*0 + 1));
			c5 = *(U_color + *(dof_el + ndof*1 + 1));
			c6 = *(U_color + *(dof_el + ndof*2 + 1));
			c7 = *(U_color + *(dof_el + ndof*3 + 1));
               	    break;
               	    case 21:
			c0 = *(U_color + *(dof_el + ndof*0 + 2));
			c1 = *(U_color + *(dof_el + ndof*1 + 2));
			c2 = *(U_color + *(dof_el + ndof*2 + 2));
			c3 = *(U_color + *(dof_el + ndof*3 + 2));
			c4 = *(U_color + *(dof_el + ndof*0 + 2));
			c5 = *(U_color + *(dof_el + ndof*1 + 2));
			c6 = *(U_color + *(dof_el + ndof*2 + 2));
			c7 = *(U_color + *(dof_el + ndof*3 + 2));
               	    break;
               	    case 22:
			c0 = *(U_color + *(dof_el + ndof*0 + 3));
			c1 = *(U_color + *(dof_el + ndof*1 + 3));
			c2 = *(U_color + *(dof_el + ndof*2 + 3));
			c3 = *(U_color + *(dof_el + ndof*3 + 3));
			c4 = *(U_color + *(dof_el + ndof*0 + 3));
			c5 = *(U_color + *(dof_el + ndof*1 + 3));
			c6 = *(U_color + *(dof_el + ndof*2 + 3));
			c7 = *(U_color + *(dof_el + ndof*3 + 3));
               	    break;
               	    case 23:
			c0 = *(U_color + *(dof_el + ndof*0 + 4));
			c1 = *(U_color + *(dof_el + ndof*1 + 4));
			c2 = *(U_color + *(dof_el + ndof*2 + 4));
			c3 = *(U_color + *(dof_el + ndof*3 + 4));
			c4 = *(U_color + *(dof_el + ndof*0 + 4));
			c5 = *(U_color + *(dof_el + ndof*1 + 4));
			c6 = *(U_color + *(dof_el + ndof*2 + 4));
			c7 = *(U_color + *(dof_el + ndof*3 + 4));
               	    break;
               	    case 30:
			c0 = 0;
			c1 = 0;
			c2 = 0;
			c3 = 0;
			c4 = 0;
			c5 = 0;
			c6 = 0;
			c7 = 0;
			if( matl_choice == matl_number )
			{
				c0 = 7;
				c1 = 7;
				c2 = 7;
				c3 = 7;
				c4 = 7;
				c5 = 7;
				c6 = 7;
				c7 = 7;
			}
               	    break;
               	    case 31:
			c0 = 0;
			c1 = 0;
			c2 = 0;
			c3 = 0;
			c4 = 0;
			c5 = 0;
			c6 = 0;
			c7 = 0;
               	    break;
               	    case 32:
			c0 = 0;
			c1 = 0;
			c2 = 0;
			c3 = 0;
			c4 = 0;
			c5 = 0;
			c6 = 0;
			c7 = 0;
			if( ele_choice == k )
			{
				c0 = 7;
				c1 = 7;
				c2 = 7;
				c3 = 7;
				c4 = 7;
				c5 = 7;
				c6 = 7;
				c7 = 7;
			}
               	    break;
        	}

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
				glColor3fv(MeshColor[c0]);
				glVertex3dv((coord_el));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord_el+9));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord_el+3));
			glEnd();
			*(norm_temp) = norm[k].face[1].x;
			*(norm_temp+1) = norm[k].face[1].y;
			*(norm_temp+2) = norm[k].face[1].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord_el+6));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord_el+3));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord_el+9));
			glEnd();
			*(norm_temp) = norm[k].face[2].x;
			*(norm_temp+1) = norm[k].face[2].y;
			*(norm_temp+2) = norm[k].face[2].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c0]);
				glVertex3dv((coord_el));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord_el+3));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord_el+12));
			glEnd();
			*(norm_temp) = norm[k].face[3].x;
			*(norm_temp+1) = norm[k].face[3].y;
			*(norm_temp+2) = norm[k].face[3].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord_el+15));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord_el+12));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord_el+3));
			glEnd();
			*(norm_temp) = norm[k].face[4].x;
			*(norm_temp+1) = norm[k].face[4].y;
			*(norm_temp+2) = norm[k].face[4].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord_el+3));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord_el+6));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord_el+15));
			glEnd();
			*(norm_temp) = norm[k].face[5].x;
			*(norm_temp+1) = norm[k].face[5].y;
			*(norm_temp+2) = norm[k].face[5].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c6]);
				glVertex3dv((coord_el+18));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord_el+15));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord_el+6));
			glEnd();
			*(norm_temp) = norm[k].face[6].x;
			*(norm_temp+1) = norm[k].face[6].y;
			*(norm_temp+2) = norm[k].face[6].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c6]);
				glVertex3dv((coord_el+18));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord_el+6));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord_el+21));
			glEnd();
			*(norm_temp) = norm[k].face[7].x;
			*(norm_temp+1) = norm[k].face[7].y;
			*(norm_temp+2) = norm[k].face[7].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord_el+9));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord_el+21));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord_el+6));
			glEnd();
			*(norm_temp) = norm[k].face[8].x;
			*(norm_temp+1) = norm[k].face[8].y;
			*(norm_temp+2) = norm[k].face[8].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord_el+21));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord_el+9));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord_el+12));
			glEnd();
			*(norm_temp) = norm[k].face[9].x;
			*(norm_temp+1) = norm[k].face[9].y;
			*(norm_temp+2) = norm[k].face[9].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c0]);
				glVertex3dv((coord_el));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord_el+12));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord_el+9));
			glEnd();
			*(norm_temp) = norm[k].face[10].x;
			*(norm_temp+1) = norm[k].face[10].y;
			*(norm_temp+2) = norm[k].face[10].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord_el+12));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord_el+15));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord_el+21));
			glEnd();
			*(norm_temp) = norm[k].face[11].x;
			*(norm_temp+1) = norm[k].face[11].y;
			*(norm_temp+2) = norm[k].face[11].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c6]);
				glVertex3dv((coord_el+18));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord_el+21));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord_el+15));
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
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord_el));
			   glVertex3dv((coord_el+3));
			   glVertex3dv((coord_el+15));
			   glVertex3dv((coord_el+12));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord_el+3));
			   glVertex3dv((coord_el+6));
			   glVertex3dv((coord_el+18));
			   glVertex3dv((coord_el+15));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord_el+21));
			   glVertex3dv((coord_el+18));
			   glVertex3dv((coord_el+6));
			   glVertex3dv((coord_el+9));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord_el+9));
			   glVertex3dv((coord_el));
			   glVertex3dv((coord_el+12));
			   glVertex3dv((coord_el+21));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord_el+12));
			   glVertex3dv((coord_el+15));
			   glVertex3dv((coord_el+18));
			   glVertex3dv((coord_el+21));
		     glEnd();

		}

		if( input_color_flag )
		{
		     c0 = 8;
		     c1 = 8;
		     c2 = 8;
		     c3 = 8;
		     c4 = 8;
		     c5 = 8;
		     c6 = 8;
		     c7 = 8;
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
				glColor3fv(MeshColor[c0]);
				glVertex3dv((coord0_el));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord0_el+9));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord0_el+3));
			glEnd();
			*(norm_temp) = norm0[k].face[1].x;
			*(norm_temp+1) = norm0[k].face[1].y;
			*(norm_temp+2) = norm0[k].face[1].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord0_el+6));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord0_el+3));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord0_el+9));
			glEnd();
			*(norm_temp) = norm0[k].face[2].x;
			*(norm_temp+1) = norm0[k].face[2].y;
			*(norm_temp+2) = norm0[k].face[2].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c0]);
				glVertex3dv((coord0_el));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord0_el+3));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord0_el+12));
			glEnd();
			*(norm_temp) = norm0[k].face[3].x;
			*(norm_temp+1) = norm0[k].face[3].y;
			*(norm_temp+2) = norm0[k].face[3].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord0_el+15));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord0_el+12));
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord0_el+3));
			glEnd();
			*(norm_temp) = norm0[k].face[4].x;
			*(norm_temp+1) = norm0[k].face[4].y;
			*(norm_temp+2) = norm0[k].face[4].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c1]);
				glVertex3dv((coord0_el+3));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord0_el+6));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord0_el+15));
			glEnd();
			*(norm_temp) = norm0[k].face[5].x;
			*(norm_temp+1) = norm0[k].face[5].y;
			*(norm_temp+2) = norm0[k].face[5].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c6]);
				glVertex3dv((coord0_el+18));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord0_el+15));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord0_el+6));
			glEnd();
			*(norm_temp) = norm0[k].face[6].x;
			*(norm_temp+1) = norm0[k].face[6].y;
			*(norm_temp+2) = norm0[k].face[6].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c6]);
				glVertex3dv((coord0_el+18));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord0_el+6));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord0_el+21));
			glEnd();
			*(norm_temp) = norm0[k].face[7].x;
			*(norm_temp+1) = norm0[k].face[7].y;
			*(norm_temp+2) = norm0[k].face[7].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord0_el+9));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord0_el+21));
				glColor3fv(MeshColor[c2]);
				glVertex3dv((coord0_el+6));
			glEnd();
			*(norm_temp) = norm0[k].face[8].x;
			*(norm_temp+1) = norm0[k].face[8].y;
			*(norm_temp+2) = norm0[k].face[8].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord0_el+21));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord0_el+9));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord0_el+12));
			glEnd();
			*(norm_temp) = norm0[k].face[9].x;
			*(norm_temp+1) = norm0[k].face[9].y;
			*(norm_temp+2) = norm0[k].face[9].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c0]);
				glVertex3dv((coord0_el));
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord0_el+12));
				glColor3fv(MeshColor[c3]);
				glVertex3dv((coord0_el+9));
			glEnd();
			*(norm_temp) = norm0[k].face[10].x;
			*(norm_temp+1) = norm0[k].face[10].y;
			*(norm_temp+2) = norm0[k].face[10].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c4]);
				glVertex3dv((coord0_el+12));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord0_el+15));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord0_el+21));
			glEnd();
			*(norm_temp) = norm0[k].face[11].x;
			*(norm_temp+1) = norm0[k].face[11].y;
			*(norm_temp+2) = norm0[k].face[11].z;
			glBegin(GL_TRIANGLES);
				glNormal3fv(norm_temp);
				glColor3fv(MeshColor[c6]);
				glVertex3dv((coord0_el+18));
				glColor3fv(MeshColor[c7]);
				glVertex3dv((coord0_el+21));
				glColor3fv(MeshColor[c5]);
				glVertex3dv((coord0_el+15));
			glEnd();
		   }
   
/* Draw the wire frame around the mesh */

		   glColor3fv( wire_color );
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord0_el+9));
			   glVertex3dv((coord0_el+6));
			   glVertex3dv((coord0_el+3));
			   glVertex3dv((coord0_el));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord0_el));
			   glVertex3dv((coord0_el+3));
			   glVertex3dv((coord0_el+15));
			   glVertex3dv((coord0_el+12));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord0_el+3));
			   glVertex3dv((coord0_el+6));
			   glVertex3dv((coord0_el+18));
			   glVertex3dv((coord0_el+15));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord0_el+21));
			   glVertex3dv((coord0_el+18));
			   glVertex3dv((coord0_el+6));
			   glVertex3dv((coord0_el+9));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord0_el+9));
			   glVertex3dv((coord0_el));
			   glVertex3dv((coord0_el+12));
			   glVertex3dv((coord0_el+21));
		     glEnd();
		     glBegin(GL_LINE_LOOP);
			   glVertex3dv((coord0_el+12));
			   glVertex3dv((coord0_el+15));
			   glVertex3dv((coord0_el+18));
			   glVertex3dv((coord0_el+21));
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
	    	fpointz = *(coord+nsd*node_number+2);
              	glBegin(GL_POINTS);
        	    glColor3fv(yellow);
               	    glVertex3f(fpointx, fpointy, fpointz);
    	   	glEnd();
	    }
    	    if( Before_gr_flag )
	    {
	    	fpointx = *(coord0+nsd*node_number);
	    	fpointy = *(coord0+nsd*node_number+1);
	    	fpointz = *(coord0+nsd*node_number+2);
              	glBegin(GL_POINTS);
        	    glColor3fv(yellow);
               	    glVertex3f(fpointx, fpointy, fpointz);
    	   	glEnd();
	    }
	}
	/*return 1;*/
}

