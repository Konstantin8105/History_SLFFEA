/*
    This program draws the drag down menus.  It works with a brick FEM code.
  
   			Last Update 6/9/01

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
#include <string.h>
#include "../brick2/br2struct.h"
#include "../brick_gr/brstrcgr.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

extern int nmat, numnp, numel, dof;
extern double *coord, *coord0;
extern double *U;
extern int *connecter;
extern BOUND bc;
extern double *force;
extern SDIM *stress_node;
extern SDIM *strain_node;
extern double *heat_el, *heat_node, *T, *Q;
extern XYZF *force_vec, *force_vec0;
extern ISTRESS *stress_color;
extern ISTRAIN *strain_color;
extern int *U_color, *T_color, *Q_color;

extern int input_flag, post_flag, color_choice,
    choice, matl_choice, node_choice, ele_choice;
extern int input_color_flag;
extern int Solid_flag, Perspective_flag, Render_flag,
    AppliedDisp_flag, AppliedForce_flag,
    Material_flag, Node_flag, Element_flag, Axes_flag;
extern int Before_flag, After_flag,
    Both_flag, Amplify_flag;
extern double amplify_factor, amplify_step, amplify_step0;

int br2set(BOUND , int *, double *, XYZF *, double *, int *, SDIM *,
        ISTRAIN *, SDIM *, ISTRESS *, double *, int *, double *, int *);

void br2ReGetparameter( void);

void br2MenuSelect(int value)
{
	int check;

	switch (value) {
	case 1:
	    color_choice = 31;
	    input_color_flag = 0;
	    AppliedForce_flag = 0;
	    AppliedDisp_flag = 0;
	    Element_flag = 0;
	    Material_flag = 0;
	    Node_flag = 1;

	    printf("\n What is the desired node number?\n");
	    scanf("%d", &node_choice);
	    if ( node_choice > numnp - 1 )
	    {
		node_choice = 0;
	    }
	    break;
	case 2:
	    color_choice = 32;
	    input_color_flag = 0;
	    AppliedForce_flag = 0;
	    AppliedDisp_flag = 0;
	    Element_flag = 1;
	    Material_flag = 0;
	    Node_flag = 0;

	    printf("\n What is the desired element number?\n");
	    scanf("%d", &ele_choice);
	    if ( ele_choice > numel - 1 )
	    {
		ele_choice = 0;
	    }
	    Solid_flag = 1;
	    break;
	case 3:
	    color_choice = 30;
	    input_color_flag = 0;
	    AppliedForce_flag = 0;
	    AppliedDisp_flag = 0;
	    Element_flag = 0;
	    Material_flag = 1;
	    Node_flag = 0;

	    printf("\n What is the desired material number?\n");
	    scanf("%d", &matl_choice);
	    if ( matl_choice > nmat - 1 )
	    {
		matl_choice = 0;
	    }
	    Solid_flag = 1;
	    break;
	case 4:
	    br2ReGetparameter();
            check = br2set( bc, connecter, force, force_vec0, Q, Q_color,
		strain_node, strain_color, stress_node, stress_color, T,
		T_color, U, U_color);
            if(!check) printf( " Problems with br2set \n");
	    break;
	case 5:
	    exit(0);
	    break;
	}

	glutPostRedisplay();

}

void br2Menu(void)
{
	glutCreateMenu(br2MenuSelect);
	glutAddMenuEntry("Jump node", 1);
	glutAddMenuEntry("Jump ele", 2);
	glutAddMenuEntry("Jump Matl", 3);
	glutAddMenuEntry("Re-Param", 4);
	glutAddMenuEntry("Quit", 5);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}
