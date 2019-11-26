/*
   This program draws the displacement and force vectors
   for quad elements. 
 
                                 San Le
 
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
#include <math.h>
#include "../quad/qdconst.h"
#include "../quad/qdstruct.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

extern GLfloat yellow[3], orange[3], orangeRed[3], red[3], green[3], 
 	violetRed[3], magenta[3], purple[3], blue[3],
	white[3], grey[3], black[3], brown[3];
extern double AxisLength_max;

void qddisp_vectors0(int displaylistnum, BOUND bc, double *coord0)
{

/* draws displacement vector on undefromed configuration */

  int i,j, dum;
  double fdum, fpointx, fpointy, fpointz;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
  fdum = AxisLength_max;
  glNewList(displaylistnum, GL_COMPILE);
  glPushAttrib(GL_LIGHTING_BIT);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, axes_ambuse);
    glLineWidth (4.0);
    glBegin(GL_LINES);
/* Draw the x displacement vectors */
	glColor3fv(white);
	for( i = 0; i < bc.num_fix[0].x; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].x);
		fpointy = *(coord0+nsd*bc.fix[i].x + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum, fpointy, fpointz); 
	}
/* Draw the y displacement vectors */
	glColor3fv(grey);
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].y);
		fpointy = *(coord0+nsd*bc.fix[i].y + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx, fpointy - fdum, fpointz); 
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(blue);
	for( i = 0; i < bc.num_fix[0].x; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].x);
		fpointy = *(coord0+nsd*bc.fix[i].x + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].y);
		fpointy = *(coord0+nsd*bc.fix[i].y + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
  glPopAttrib();
  glEndList();
}

void qdforce_vectors0(int displaylistnum, BOUND bc, double *coord0,
	XYF *force_vec )
{

/* draws force vector on undefromed configuration */

  int i,j, dum;
  double fpointx, fpointy, fpointz, fx, fy, fz;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
  glNewList(displaylistnum, GL_COMPILE);
  glPushAttrib(GL_LIGHTING_BIT);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, axes_ambuse);
    glLineWidth (4.0);
    glBegin(GL_LINES);
/* Draw the force vectors axis */
	glColor3fv(white);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
		fx = force_vec[i].x; fy = force_vec[i].y;
			fz = 0.0;
		fpointx = *(coord0+nsd*bc.force[i]);
		fpointy = *(coord0+nsd*bc.force[i] + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fx, fy, fz); 
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(red);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
		fpointx = *(coord0+nsd*bc.force[i]);
		fpointy = *(coord0+nsd*bc.force[i] + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
  glPopAttrib();
  glEndList();
}

void qddisp_vectors(BOUND bc, double *coord)
{

/* draws displacement vector on deformed configuration */

  int i,j, dum;
  double fdum, fpointx, fpointy, fpointz;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
  fdum = AxisLength_max;
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, axes_ambuse);
    glLineWidth (4.0);
    glBegin(GL_LINES);
/* Draw the x displacement vectors */
	glColor3fv(white);
	for( i = 0; i < bc.num_fix[0].x; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].x);
		fpointy = *(coord+nsd*bc.fix[i].x + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum, fpointy, fpointz); 
	}
/* Draw the y displacement vectors */
	glColor3fv(grey);
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].y);
		fpointy = *(coord+nsd*bc.fix[i].y + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx, fpointy - fdum, fpointz); 
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(blue);
	for( i = 0; i < bc.num_fix[0].x; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].x);
		fpointy = *(coord+nsd*bc.fix[i].x + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].y);
		fpointy = *(coord+nsd*bc.fix[i].y + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
}

void qdforce_vectors(BOUND bc, double *coord, XYF *force_vec )
{

/* draws force vector on deformed configuration */

  int i,j, dum;
  double fpointx, fpointy, fpointz, fx, fy, fz;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, axes_ambuse);
    glLineWidth (4.0);
    glBegin(GL_LINES);
/* Draw the force vectors axis */
	glColor3fv(white);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
		fx = force_vec[i].x; fy = force_vec[i].y;
			fz = 0.0;
		fpointx = *(coord+nsd*bc.force[i]);
		fpointy = *(coord+nsd*bc.force[i] + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fx, fy, fz); 
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(red);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
		fpointx = *(coord+nsd*bc.force[i]);
		fpointy = *(coord+nsd*bc.force[i] + 1);
		fpointz = 0.0;
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
}

