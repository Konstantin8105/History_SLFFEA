/*
    This program contains the mesh display routine for the FEM GUI
    for plate elements.
  
   			Last Update 6/30/99

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
 */
#include <stdio.h>
#include <stdlib.h>
#include "../plate/plconst.h"
#include "../plate/plstruct.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void plmeshdraw(void);
void plrender(void);
void pldisp_vectors(BOUND , double *, double *);
void plforce_vectors(BOUND , double *, double *, ZPhiF *);

extern double *coord, *zcoord;
extern BOUND bc;

extern double left_right, up_down, in_out, left_right0,
	up_down0, in_out0, xAngle, yAngle, zAngle;
extern GLuint AxesList, DispList, ForceList;   /* Display lists */
extern ZPhiF *force_vec;
extern int Render_flag, AppliedDisp_flag, AppliedForce_flag,
    Axes_flag, Before_flag, After_flag; 

void plMeshDisplay(void)
{
    	glClear (GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);

    	glLoadIdentity ();  /*  clear the matrix    */

        glTranslatef (left_right, up_down, in_out);
        glRotatef (xAngle, 1, 0, 0);
        glRotatef (yAngle, 0, 1, 0);
        glRotatef (zAngle, 0, 0, 1);

    	glPointSize(8);
	if(Axes_flag)
		glCallList(AxesList);
	if(AppliedDisp_flag)
	{
		if(Before_flag )
			glCallList(DispList);
		if(After_flag )
  			pldisp_vectors(bc, coord, zcoord);
	}
	if(AppliedForce_flag)
	{
		if(Before_flag )
			glCallList(ForceList);
		if(After_flag )
  			plforce_vectors(bc, coord, zcoord, force_vec);
	}
    	glPushMatrix ();
	glLineWidth (2.0);
	if(Render_flag)
		plrender();
    	else
		plmeshdraw();
    	glPopMatrix ();
	glFlush();
  	glutSwapBuffers();
}
