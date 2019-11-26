/*
    This program contains the control display routine for the FEM GUI
    for beam elements.
  
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../beam/bmconst.h"
#include "../beam/bmstruct.h"
#include "bmstrcgr.h"
#include "../../common_gr/control.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

extern GLfloat MeshColor[boxnumber+5][3];

extern int ControlDiv_y[rowdim + 2], ControlDiv_x[rowdim + 2];
extern char ControlText[][90];
extern int del_height;
extern int del_width;
extern double ratio, ratio2;
extern int boxMove_x, boxMove_y, boxTextMove_x, textMove_x, textMove_y[rowdim+2];

extern int current_width, current_height;
extern int control_height, control_width, mesh_height, mesh_width;
extern int stress_flag, strain_flag, stress_strain, disp_flag, angle_flag;

extern GLfloat yellow[3], orange[3], orangeRed[3], red[3], green[3],
        violetRed[3], magenta[3], purple[3], blue[3],
        white[3], grey[3], black[3];

extern char RotateData[3][10];
extern char MoveData[3][10];
extern char AmplifyData[10];
extern char BoxData[2*boxnumber+2][14];
extern char BoxText[10];

extern int Color_flag[rowdim];

void printText(const char *);

void bmControlDisplay(void)
{
	int i, j, dum = 0, dum2 = 0;
    	va_list args;
    	GLfloat font_scale = 119.05 + 33.33;

    	glClear(GL_COLOR_BUFFER_BIT);
    	glMatrixMode(GL_MODELVIEW);
    	glPushMatrix();
    	glColor3fv(white);

/* Draw the Text */

/* This draws the View Option Text */

/* Calculate the translation in Y for the text */

    	for( i = 0 ; i < rowdim + 2 ; ++i)
	{
		textMove_y[i] = current_height -
			2*textHeight - (int)(ratio2*textHeightDiv*dum);
		++dum;
	}

    	for( i = 0 ; i < 2 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    	   	printText( ControlText[2*i] );
    	   	printText( ControlText[2*i+1] );
	}

/* Change displayed text for Node ID, Element ID, and Material */

	if( Color_flag[3])
	{
		strcpy(ControlText[6],"Next   Last  ");
	}
	else
	{
		strcpy(ControlText[6]," On          ");
	}
	if( Color_flag[4])
	{
		strcpy(ControlText[8],"Next   Last  ");
	}
	else
	{
		strcpy(ControlText[8]," On          ");
	}
	if( Color_flag[5])
	{
		strcpy(ControlText[10],"Next   Last  ");
	}
	else
	{
		strcpy(ControlText[10]," On          ");
	}

    	for( i = 2 ; i < 9 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    		glColor3fv(white);
		if( Color_flag[i])
		{
    			glColor3fv(yellow);
		}
    	   	printText( ControlText[2*i] );
    		glColor3fv(white);
    	   	printText( ControlText[2*i+1] );
	}

/* This draws the Rotation and Move Text */

    	glColor3fv(white);
    	for( i = 9 ; i < 11 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    	   	printText( ControlText[2*i] );
    	   	printText( ControlText[2*i+1] );
	}
/* This draws the Rotation  */
    	for( i = 11 ; i < 14 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    	   	printText( ControlText[2*i] );
    	   	printText( ControlText[2*i+1] );
    	   	printText( RotateData[dum2]);
		++dum2;
	}
	dum2 = 0;
        for( i = 14 ; i < 16 ; ++i)
        {
                glLoadIdentity();
                glTranslatef (0,textMove_y[i],0);
                glScalef( ratio, ratio, 1.0);
                printText( ControlText[2*i] );
                printText( ControlText[2*i+1] );
        }
/* This draws the Move Text */
    	for( i = 16 ; i < 19 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    	   	printText( ControlText[2*i] );
    	   	printText( ControlText[2*i+1] );
    	   	printText( MoveData[dum2]);
		++dum2;
	}
    	for( i = 19 ; i < 23 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    	   	printText( ControlText[2*i] );
    	   	printText( ControlText[2*i+1] );
	}

/* This draws the Deformation Text */
    	glColor3fv(white);
    	for( i = 23 ; i < 26 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    		glColor3fv(white);
		if( Color_flag[i])
		{
    			glColor3fv(yellow);
		}
    	   	printText( ControlText[2*i] );
    		glColor3fv(white);
    	   	printText( ControlText[2*i+1] );
	}
/* This draws the Amplification Data  */

       	glLoadIdentity();
       	glTranslatef (0,textMove_y[26],0);
       	glScalef( ratio, ratio, 1.0);
       	printText( ControlText[52] );
       	printText( ControlText[53] );
       	printText( AmplifyData );
	++dum2;
	
/* This draws the Engineering Analysis Option Text */
/* Stress and Strain */
    	glColor3fv(white);
        for( i = 27 ; i < 29 ; ++i)
        {
                glLoadIdentity();
                glTranslatef (0,textMove_y[i],0);
                glScalef( ratio, ratio, 1.0);
                printText( ControlText[2*i] );
                printText( ControlText[2*i+1] );
        }
    	for( i = 29 ; i < 34 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    		glColor3fv(white);
		if( Color_flag[i])
		{
			if( stress_flag)
		  	{
    		  		glColor3fv(yellow);
		  	}
		}
    	   	printText( ControlText[2*i] );
    		glColor3fv(white);
		if( Color_flag[i])
		{
		  	if( strain_flag)
		  	{
    		  		glColor3fv(yellow);
		  	}
		}
    	   	printText( ControlText[2*i+1] );
	}

/* Displacement */
    	glColor3fv(white);
        for( i = 34 ; i < 36 ; ++i)
        {
                glLoadIdentity();
                glTranslatef (0,textMove_y[i],0);
                glScalef( ratio, ratio, 1.0);
                printText( ControlText[2*i] );
                printText( ControlText[2*i+1] );
        }
    	for( i = 36 ; i < 39 ; ++i)
	{
    	   	glLoadIdentity();
    	   	glTranslatef (0,textMove_y[i],0);
    	   	glScalef( ratio, ratio, 1.0);
    		glColor3fv(white);
		if( Color_flag[i])
		{
			if( disp_flag)
		  	{
    		  		glColor3fv(yellow);
		  	}
		}
    	   	printText( ControlText[2*i] );
    		glColor3fv(white);
		if( Color_flag[i])
		{
		  	if( angle_flag)
		  	{
    		  		glColor3fv(yellow);
		  	}
		}
    	   	printText( ControlText[2*i+1] );
	}

/* Text label for Color Scale Boxes */

    	glColor3fv(white);
	boxTextMove_x = (int)(ratio2*boxTextMove_x0);
        for( i = 27 ; i < 28 ; ++i)
        {
                glLoadIdentity();
                glTranslatef (boxTextMove_x,textMove_y[i],0);
                glScalef( ratio, ratio, 1.0);
                printText( BoxText );
        }

/* Text for Color Scale Boxes */

        dum = 0;
    	glColor3fv(white);
        for( i = 29 ; i < rowdim + 2 ; i += 1)
        {
                glLoadIdentity();
                glTranslatef (boxTextMove_x,textMove_y[i] + textHeight,0);
                glScalef( ratio, ratio, 1.0);
                printText( BoxData[dum] );
                ++dum;
        }

/* Begin Drawing the Color Scale Boxes */
	
	del_width = current_width - ratio2*control_width0;
	del_height = current_height - ratio2*control_height0;
    	boxMove_x = (int)(ratio2*(control_width0 - left_indent));
    	boxMove_y = del_height + (int)(ratio2*bottom_indent);

	glLoadIdentity();
	glTranslatef (boxMove_x,boxMove_y,0);
	glScalef( ratio2, ratio2, 1.0);
	glColor3fv(MeshColor[0]);
	glRects(0,0,boxdim,boxdim);
	glTranslatef (0,boxHeight,0);
	glColor3fv(MeshColor[1]);
	glRects(0,0,boxdim,boxdim);
	glTranslatef (0,boxHeight,0);
	glColor3fv(MeshColor[2]);
	glRects(0,0,boxdim,boxdim);
	glTranslatef (0,boxHeight,0);
	glColor3fv(MeshColor[3]);
	glRects(0,0,boxdim,boxdim);
	glTranslatef (0,boxHeight,0);
	glColor3fv(MeshColor[4]);
	glRects(0,0,boxdim,boxdim);
	glTranslatef (0,boxHeight,0);
	glColor3fv(MeshColor[5]);
	glRects(0,0,boxdim,boxdim);
	glTranslatef (0,boxHeight,0);
	glColor3fv(MeshColor[6]);
	glRects(0,0,boxdim,boxdim);
	glTranslatef (0,boxHeight,0);
	glColor3fv(MeshColor[7]);
	glRects(0,0,boxdim,boxdim);

	glutSwapBuffers();
}

