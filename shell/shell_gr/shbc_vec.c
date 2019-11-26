/*
   This program draws the displacement and force vectors
   for shell elements. 
 
                                 San Le
 
  			Last Update 5/14/00
 */

#if WINDOWS
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../shell/shconst.h"
#include "../shell/shstruct.h"
#include "../../common_gr/control.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

extern GLfloat yellow[3], orange[3], orangeRed[3], red[3], green[3], 
 	violetRed[3], magenta[3], purple[3], blue[3],
	white[3], grey[3], black[3], brown[3];
extern double AxisLength_max;

void shdisp_vectors0(int displaylistnum, BOUND bc, double *coord0)
{

/* draws displacement vector on undefromed configuration */

  int i,j, dum;
  double fdum, fdum2, fpointx, fpointy, fpointz;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
  fdum = AxisLength_max;
  fdum2 = .5*AxisLength_max;
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
		fpointz = *(coord0+nsd*bc.fix[i].x + 2);
		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum, fpointy, fpointz);
	}
/* Draw the y displacement vectors */
	glColor3fv(grey);
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].y);
		fpointy = *(coord0+nsd*bc.fix[i].y + 1);
		fpointz = *(coord0+nsd*bc.fix[i].y + 2);
		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx, fpointy - fdum, fpointz);
	}
/* Draw the z displacement vectors */
	glColor3fv(black);
	for( i = 0; i < bc.num_fix[0].z; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].z);
		fpointy = *(coord0+nsd*bc.fix[i].z + 1);
		fpointz = *(coord0+nsd*bc.fix[i].z + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx, fpointy, fpointz - fdum); 
	}
    glEnd();
    glLineWidth (2.0);
    glBegin(GL_LINES);
/* Draw the x angle vectors */
	glColor3fv(white);
	for( i = 0; i < bc.num_fix[0].phix; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].phix);
		fpointy = *(coord0+nsd*bc.fix[i].phix + 1);
		fpointz = *(coord0+nsd*bc.fix[i].phix + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy + fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy + fdum2, fpointz); 
		glVertex3f( fpointx, fpointy, fpointz); 
	}
/* Draw the y angle vectors */
	glColor3fv(grey);
	for( i = 0; i < bc.num_fix[0].phiy; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].phiy);
		fpointy = *(coord0+nsd*bc.fix[i].phiy + 1);
		fpointz = *(coord0+nsd*bc.fix[i].phiy + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx + fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx + fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx, fpointy, fpointz); 
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(blue);
	for( i = 0; i < bc.num_fix[0].x; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].x);
		fpointy = *(coord0+nsd*bc.fix[i].x + 1);
		fpointz = *(coord0+nsd*bc.fix[i].x + 2);
		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].y);
		fpointy = *(coord0+nsd*bc.fix[i].y + 1);
		fpointz = *(coord0+nsd*bc.fix[i].y + 2);
		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].z; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].z);
		fpointy = *(coord0+nsd*bc.fix[i].z + 1);
		fpointz = *(coord0+nsd*bc.fix[i].z + 2);
		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].phix; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].phix);
		fpointy = *(coord0+nsd*bc.fix[i].phix + 1);
		fpointz = *(coord0+nsd*bc.fix[i].phix + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].phiy; ++i)
	{
		fpointx = *(coord0+nsd*bc.fix[i].phiy);
		fpointy = *(coord0+nsd*bc.fix[i].phiy + 1);
		fpointz = *(coord0+nsd*bc.fix[i].phiy + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
  glPopAttrib();
  glEndList();
}

void shforce_vectors0(int displaylistnum, BOUND bc, double *coord0,
	XYZPhiF *force_vec )
{

/* draws force vector on undefromed configuration */

  int i,j, dum;
  double fpointx, fpointy, fpointz, fx, fy, fz, fdum, fdum2, fdum3;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
  glNewList(displaylistnum, GL_COMPILE);
  glPushAttrib(GL_LIGHTING_BIT);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, axes_ambuse);
    glLineWidth (4.0);
    glBegin(GL_LINES);
/* Draw the force vectors */
	glColor3fv(white);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
	    fx = force_vec[i].x; fy = force_vec[i].y;
	    	fz = force_vec[i].z;
	    fpointx = *(coord0+nsd*bc.force[i]);
	    fpointy = *(coord0+nsd*bc.force[i] + 1);
	    fpointz = *(coord0+nsd*bc.force[i] + 2);
	    fdum = fabs(fx-fpointx);
	    fdum += fabs(fy-fpointy);
	    fdum += fabs(fz-fpointz);
	    if( fdum > SMALL)
	    {
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fx, fy, fz); 
	    }
	}
    glEnd();
    glLineWidth (2.0);
    glBegin(GL_LINES);
/* Draw the moment vectors */
	glColor3fv(black);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
	    fx = force_vec[i].phix; fy = force_vec[i].phiy;
	    fpointx = *(coord0+nsd*bc.force[i]);
	    fpointy = *(coord0+nsd*bc.force[i] + 1);
	    fpointz = *(coord0+nsd*bc.force[i] + 2);
	    fdum = fx-fpointx;
	    fdum2 = fy-fpointy;
	    fdum3 = fabs(fdum) + fabs(fdum2);
	    if( fdum3 > SMALL)
	    {
		fdum *= .2;
		fdum2 *= .2;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fx + fdum2, fy - fdum, fpointz); 
		glVertex3f( fx + fdum2, fy - fdum, fpointz); 
		glVertex3f( fx - fdum2, fy + fdum, fpointz); 
		glVertex3f( fx - fdum2, fy + fdum, fpointz); 
      		glVertex3f( fpointx, fpointy, fpointz);
	    }
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(red);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
		fpointx = *(coord0+nsd*bc.force[i]);
		fpointy = *(coord0+nsd*bc.force[i] + 1);
		fpointz = *(coord0+nsd*bc.force[i] + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
  glPopAttrib();
  glEndList();
}

void shdisp_vectors(BOUND bc, double *coord)
{

/* draws displacement vector on deformed configuration */

  int i,j, dum;
  double fdum, fdum2, fpointx, fpointy, fpointz;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
  fdum = AxisLength_max;
  fdum2 = .5*AxisLength_max;
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, axes_ambuse);
    glLineWidth (4.0);
    glBegin(GL_LINES);
/* Draw the x displacement vectors */
	glColor3fv(white);
	for( i = 0; i < bc.num_fix[0].x; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].x);
		fpointy = *(coord+nsd*bc.fix[i].x + 1);
		fpointz = *(coord+nsd*bc.fix[i].x + 2);
		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum, fpointy, fpointz);
	}
/* Draw the y displacement vectors */
	glColor3fv(grey);
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].y);
		fpointy = *(coord+nsd*bc.fix[i].y + 1);
		fpointz = *(coord+nsd*bc.fix[i].y + 2);
		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx, fpointy - fdum, fpointz);
	}
/* Draw the z displacement vectors */
	glColor3fv(black);
	for( i = 0; i < bc.num_fix[0].z; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].z);
		fpointy = *(coord+nsd*bc.fix[i].z + 1);
		fpointz = *(coord+nsd*bc.fix[i].z + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx, fpointy, fpointz - fdum); 
	}
    glEnd();
    glLineWidth (2.0);
    glBegin(GL_LINES);
/* Draw the x angle vectors */
	glColor3fv(white);
	for( i = 0; i < bc.num_fix[0].phix; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].phix);
		fpointy = *(coord+nsd*bc.fix[i].phix + 1);
		fpointz = *(coord+nsd*bc.fix[i].phix + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy + fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy + fdum2, fpointz); 
		glVertex3f( fpointx, fpointy, fpointz); 
	}
/* Draw the y angle vectors */
	glColor3fv(grey);
	for( i = 0; i < bc.num_fix[0].phiy; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].phiy);
		fpointy = *(coord+nsd*bc.fix[i].phiy + 1);
		fpointz = *(coord+nsd*bc.fix[i].phiy + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx - fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx + fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx + fdum2, fpointy - fdum2, fpointz); 
		glVertex3f( fpointx, fpointy, fpointz); 
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(blue);
	for( i = 0; i < bc.num_fix[0].x; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].x);
		fpointy = *(coord+nsd*bc.fix[i].x + 1);
		fpointz = *(coord+nsd*bc.fix[i].x + 2);
		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].y; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].y);
		fpointy = *(coord+nsd*bc.fix[i].y + 1);
		fpointz = *(coord+nsd*bc.fix[i].y + 2);
		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].z; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].z);
		fpointy = *(coord+nsd*bc.fix[i].z + 1);
		fpointz = *(coord+nsd*bc.fix[i].z + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].phix; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].phix);
		fpointy = *(coord+nsd*bc.fix[i].phix + 1);
		fpointz = *(coord+nsd*bc.fix[i].phix + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
	}
	for( i = 0; i < bc.num_fix[0].phiy; ++i)
	{
		fpointx = *(coord+nsd*bc.fix[i].phiy);
		fpointy = *(coord+nsd*bc.fix[i].phiy + 1);
		fpointz = *(coord+nsd*bc.fix[i].phiy + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
}

void shforce_vectors(BOUND bc, double *coord, XYZPhiF *force_vec )
{

/* draws force vector on deformed configuration */

  int i,j, dum;
  double fpointx, fpointy, fpointz, fx, fy, fz, fdum, fdum2, fdum3;
  GLfloat axes_ambuse[] =   { 0.5, 0.0, 0.0, 1.0 };
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, axes_ambuse);
    glLineWidth (4.0);
    glBegin(GL_LINES);
/* Draw the force vectors */
	glColor3fv(white);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
	    fx = force_vec[i].x; fy = force_vec[i].y;
	    	fz = force_vec[i].z;
	    fpointx = *(coord+nsd*bc.force[i]);
	    fpointy = *(coord+nsd*bc.force[i] + 1);
	    fpointz = *(coord+nsd*bc.force[i] + 2);
	    fdum = fabs(fx-fpointx);
	    fdum += fabs(fy-fpointy);
	    fdum += fabs(fz-fpointz);
	    if( fdum > SMALL)
	    {
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fx, fy, fz); 
	    }
	}
    glEnd();
    glLineWidth (2.0);
    glBegin(GL_LINES);
/* Draw the moment vectors */
	glColor3fv(black);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
	    fx = force_vec[i].phix; fy = force_vec[i].phiy;
	    fpointx = *(coord+nsd*bc.force[i]);
	    fpointy = *(coord+nsd*bc.force[i] + 1);
	    fpointz = *(coord+nsd*bc.force[i] + 2);
	    fdum = fx-fpointx;
	    fdum2 = fy-fpointy;
	    fdum3 = fabs(fdum) + fabs(fdum2);
	    if( fdum3 > SMALL)
	    {
		fdum *= .2;
		fdum2 *= .2;
      		glVertex3f( fpointx, fpointy, fpointz);
		glVertex3f( fx + fdum2, fy - fdum, fpointz); 
		glVertex3f( fx + fdum2, fy - fdum, fpointz); 
		glVertex3f( fx - fdum2, fy + fdum, fpointz); 
		glVertex3f( fx - fdum2, fy + fdum, fpointz); 
      		glVertex3f( fpointx, fpointy, fpointz);
	    }
	}
    glEnd();
    glPointSize(8);
    glBegin(GL_POINTS);
	glColor3fv(red);
	for( i = 0; i < bc.num_force[0]; ++i)
	{
		fpointx = *(coord+nsd*bc.force[i]);
		fpointy = *(coord+nsd*bc.force[i] + 1);
		fpointz = *(coord+nsd*bc.force[i] + 2);
      		glVertex3f( fpointx, fpointy, fpointz);
	}
    glEnd();
}

