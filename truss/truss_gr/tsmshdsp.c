/*
    This program contains the mesh display routine for the FEM GUI
    for truss elements.
  
   			Last Update 5/14/00

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
#include "../truss/tsconst.h"
#include "../truss/tsstruct.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void tsmeshdraw(void);
void tsdisp_vectors(BOUND , double *);
void tsforce_vectors(BOUND , double *, XYZF *);

extern double *coord;
extern BOUND bc;

extern double left_right, up_down, in_out, left_right0,
	up_down0, in_out0, xAngle, yAngle, zAngle;
extern GLuint AxesList, DispList, ForceList;   /* Display lists */
extern XYZF *force_vec;
extern int Render_flag, AppliedDisp_flag, AppliedForce_flag,
    Axes_flag, Before_flag, After_flag; 

void tsMeshDisplay(void)
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
  			tsdisp_vectors(bc, coord);
	}
	if(AppliedForce_flag)
	{
		if(Before_flag )
			glCallList(ForceList);
		if(After_flag )
  			tsforce_vectors(bc, coord, force_vec);
	}
    	glPushMatrix ();
	glLineWidth (2.0);
	tsmeshdraw();
    	glPopMatrix ();
	glFlush();
  	glutSwapBuffers();
}

