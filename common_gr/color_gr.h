/*
    This is the include file "color_gr.h" which contains the colors
    used in the graphics programs for the FEM code. 

                Updated 4/22/00

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "control.h"

/* glut header files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#define SLFFEA                0    /* For Colors based on SLFFEA */
#define ANSYS                 1    /* For Colors based on ANSYS */
#define ANSYS_2               0    /* For Colors based on ANSYS_2*/

#if SLFFEA
GLfloat MeshColor[boxnumber+5][3] = {
	0.0, 0.0, 1.0,
        0.627, 0.1255, 0.941,
        1.0, 0.0, 1.0,
        0.8157, 0.1255, 0.5647,
        1.0, 0, 0,
        1.0, 0.55, 0,
        1.0, 0.647, 0,
        1.0, 1.0, 0,
        0.0, 1.0, 0.0,
 	1.0, 1.0, 1.0,
	0.75, 0.75, 0.75,
 	0.5, 0.5, 0.5,
	0.0, 0.0, 0.0
};
#endif

#if ANSYS
GLfloat MeshColor[boxnumber+5][3] = {
        0.0,  0.0,  1.0,
        0.4, 0.8,  1.0,
        0.0,  1.0,  0.85,
        0.0,  1.0,  0.0,
        1.0,  1.0,  0,
        1.0,  0.70,  0,
        1.0,  0.54,  0,
        1.0, 0, 0,
        1.0, 0.0, 1.0,
 	1.0, 1.0, 1.0,
	0.75, 0.75, 0.75,
 	0.5, 0.5, 0.5,
	0.0, 0.0, 0.0
};
#endif

#if ANSYS_2
GLfloat MeshColor[boxnumber+5][3] = {
        0.0,  0.0,  1.0,
        0.3, 0.75,  0.9,
        0.0,  1.0,  0.8,
        0.0,  1.0,  0.0,
        0.7,  1.0,  0,
        1.0,  1.0,  0,
        1.0,  0.7,  0,
        1.0, 0, 0,
        1.0, 0.0, 1.0,
 	1.0, 1.0, 1.0,
	0.75, 0.75, 0.75,
 	0.5, 0.5, 0.5,
	0.0, 0.0, 0.0
};
#endif


GLfloat RenderColor[4] = { 0.5, 1.0, 1.0, 0.0};

GLfloat white[3] = {  1.0, 1.0, 1.0 };
GLfloat grey[3] = { 0.75, 0.75, 0.75 };
GLfloat darkGrey[3] = { 0.5, 0.5, 0.5 };
GLfloat black[3] = { 0.0, 0.0, 0.0 };
GLfloat green[3] = { 0.0, 1.0, 0.0 };
GLfloat brown[3] = { 1.0, 0.2545, 0.2545 };
GLfloat wire_color[3] = { 0.0, 0.0, 0.0 };

GLfloat yellow[3] = { 1.0, 1.0, 0 };
GLfloat orange[3] = { 1.0, 0.647, 0 };
GLfloat orangeRed[3] = { 1.0, 0.55, 0 };
GLfloat red[3] = { 1.0, 0, 0 }; 
GLfloat violetRed[3] = { 0.8157, 0.1255, 0.5647 };
GLfloat magenta[3] = { 1.0, 0.0, 1.0 };
GLfloat purple[3] = { 0.627, 0.1255, 0.941 };
GLfloat blue[3] = { 0.0, 0.0, 1.0 };

GLfloat yellowRed[3] = { 1.0,  0.70,  0 };
GLfloat greenYellow[3] = { 0.7, 1.0, 0.0 };
GLfloat blueGreen[3] = { 0.0,  1.0,  0.85 };
GLfloat aqua[3] = { 0.433, 0.837,  1.0 };
