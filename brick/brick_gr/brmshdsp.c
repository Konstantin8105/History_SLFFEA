/*
    This program contains the mesh display routine for the FEM GUI
    for brick elements.
  
   			Last Update 5/27/01

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
 */

#if WINDOWS
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "../brick/brconst.h"
#if BRICK1
#include "../brick/brstruct.h"
#endif
#if BRICK2
#include "../brick2/br2struct.h"
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void brmeshdraw(void);
void brrender(void);
void brdisp_vectors(BOUND , double *);
void brforce_vectors(BOUND , double *, XYZF *);

extern double *coord;
extern BOUND bc;

extern double left_right, up_down, in_out, left_right0,
	up_down0, in_out0, xAngle, yAngle, zAngle;
extern GLuint AxesList, DispList, ForceList;   /* Display lists */
extern XYZF *force_vec;
extern int Render_flag, AppliedDisp_flag, AppliedForce_flag,
    Axes_flag, Before_flag, After_flag; 
extern int CrossSection_flag;

void AxesNumbers2(void);

void AxesNumbers(void);

void AxesLabel(void);

void CrossSetionPlaneDraw(void);

void brMeshDisplay(void)
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
  			brdisp_vectors(bc, coord);
	}
	if(AppliedForce_flag)
	{
		if(Before_flag )
			glCallList(ForceList);
		if(After_flag )
  			brforce_vectors(bc, coord, force_vec);
	}
	if(CrossSection_flag)
	{
		CrossSetionPlaneDraw();
	}
    	glPushMatrix ();
	glLineWidth (2.0);
	if(Render_flag)
		brrender();
    	else
		brmeshdraw();
	if(Axes_flag)
	{
		/*AxesNumbers();*/
		AxesNumbers2();
	}
    	glPopMatrix ();
	glFlush();
  	glutSwapBuffers();
}

