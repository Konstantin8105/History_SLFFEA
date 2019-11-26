/*
    This program contains the control mouse routine for the FEM GUI
    for shell elements.
  
                        Last Update 10/15/06

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
#include "../shell/shconst.h"
#include "../shell/shstruct.h"
#include "shstrcgr.h"
#include "../../common_gr/control.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/****** FEA globals ******/

extern int dof, nmat, numnp, numel;
extern double *coord, *coord0;
extern double *U, *Uz_fib;
extern int *connecter;
extern BOUND bc;
extern XYZPhiF *force_vec, *force_vec0;

/* Global variables for the mesh color and nodal data */

extern int *el_matl_color;
extern MATL *matl_crtl;

/****** EXTERNAL VARIABLES ********/

extern NORM *norm, *norm0;

extern int ControlDiv_y[rowdim + 2], ControlDiv_x[rowdim + 2];
extern int boxMove_x, boxMove_y, boxTextMove_x, textMove_x, textMove_y[rowdim+2];
extern int textDiv_xa, textDiv_xb;
extern double matl_choicef, node_choicef, ele_choicef;

extern int ControlWindow, MeshWindow;
extern double step_sizex, step_sizey, step_sizez;
extern int control_height, control_width, mesh_height, mesh_width;
extern double xAngle, yAngle, zAngle;
extern double amplify_factor, amplify_step, amplify_step0;
extern double left_right, up_down, in_out;
extern double left_right0, up_down0, in_out0;
extern double ortho_left, ortho_right, ortho_top, ortho_bottom,
        ortho_left0, ortho_right0, ortho_top0, ortho_bottom0;
extern int ortho_redraw_flag;

extern int input_flag, post_flag, color_choice, matl_choice, node_choice, ele_choice;
extern int input_color_flag;
extern int Solid_flag, Perspective_flag, Render_flag, AppliedDisp_flag,
	AppliedForce_flag, Material_flag, Node_flag, Element_flag, Axes_flag;
extern int Before_flag, After_flag, Both_flag, Amplify_flag; 
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
extern double Ux_div[boxnumber+1], Uy_div[boxnumber+1], Uz_div[boxnumber+1],
	Uphi_x_div[boxnumber+1], Uphi_y_div[boxnumber+1];
extern SDIM stress_div[boxnumber+1], strain_div[boxnumber+1];

int shnormal_vectors (int *, double *, NORM * );

void ScreenShot( int , int );

void shControlMouse(int button, int state, int x, int y)
{
	int i, check, dum1, dum2;
	double fpointx, fpointy, fpointz;
  	if (button == GLUT_LEFT_BUTTON)
  	{
		if ( x < textDiv_xa )
     		{

/* These are for the View Option Keys */

     			if ( y >= ControlDiv_y[3] && y < ControlDiv_y[4] )
     			{
/* Solid Turned On */
				Solid_flag = 1;
     			}
     			if ( y >= ControlDiv_y[4] && y < ControlDiv_y[5] )
     			{
/* Node ID Turned On increment up */
				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 1;
                		node_choicef += .5;
                		if ( node_choicef > 1.0 )
                		{
                			node_choicef = 0.0;
                		}
                		node_choice += (int)node_choicef;
                		if ( node_choice > 2*numnp - 1 )
                		{
                        		node_choice = 0;
                		}
				color_choice = 31;
				Render_flag = 0;
				strain_flag = 0;
				stress_flag = 0;
     				if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}
     			}
     			if ( y >= ControlDiv_y[5] && y < ControlDiv_y[6] )
     			{
/* Element ID Turned On increment up */
				/*ScreenShot( control_width, control_height );*/

				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Element_flag = 1;
				Material_flag = 0;
				Node_flag = 0;
                		ele_choicef += .5;
                		if ( ele_choicef > 1.0 )
                		{
                			ele_choicef = 0.0;
                		}
                		ele_choice += (int)ele_choicef;
                		if ( ele_choice > numel - 1 )
                		{
                        		ele_choice = 0;
                		}
				color_choice = 32;
				Render_flag = 0;
				Solid_flag = 1;
				strain_flag = 0;
				stress_flag = 0;
     				if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}
     			}
     			if ( y >= ControlDiv_y[6] && y < ControlDiv_y[7] )
     			{
/* Material Turned On increment up */
				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Element_flag = 0;
				Node_flag = 0;
				Material_flag = 1;
                		matl_choicef += .5;
                		if ( matl_choicef > 1.0 )
                		{
                			matl_choicef = 0.0;
                		}
                		matl_choice += (int)matl_choicef;
                		if ( matl_choice > nmat - 1 )
                		{
                        		matl_choice = 0;
                		}
				color_choice = 30;
				Render_flag = 0;
				Solid_flag = 1;
				strain_flag = 0;
				stress_flag = 0;
     				if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}
     			}
     			if ( y >= ControlDiv_y[7] && y < ControlDiv_y[8] )
     			{
/* Fixed Disp Turned On */
				/*AppliedForce_flag = 0;*/
				AppliedDisp_flag = 1;
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 0;
				Render_flag = 0;
     				/*if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}*/
     			}
     			if ( y >= ControlDiv_y[8] && y < ControlDiv_y[9] )
     			{
/* Applied Force Turned On */
				AppliedForce_flag = 1;
				/*AppliedDisp_flag = 0;*/
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 0;
				Render_flag = 0;
     				/*if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}*/
     			}
     			if ( y >= ControlDiv_y[9] && y < ControlDiv_y[10] )
     			{
/* Axes Turned On */
				Axes_flag = 1;
     			}
/* These are for the Rotation Keys */
     			if ( y >= ControlDiv_y[12] && y < ControlDiv_y[13] )
     			{
/* Rotate -x */
				xAngle -= 5.0;
     			}
     			if ( y >= ControlDiv_y[13] && y < ControlDiv_y[14] )
     			{
/* Rotate -y */
				yAngle -= 5.0;
     			}
     			if ( y >= ControlDiv_y[14] && y < ControlDiv_y[15] )
     			{
/* Rotate -z */
				zAngle -= 5.0;
     			}
     			if ( y >= ControlDiv_y[15] && y < ControlDiv_y[16] )
     			{
/* Reset Angles */
				xAngle = 0.0; yAngle = 0.0; zAngle = 0.0;
     			}
/* These are for the Move Keys */
     			if ( y >= ControlDiv_y[17] && y < ControlDiv_y[18] )
     			{
/* Move -x */
				left_right -= step_sizex;
     			}
     			if ( y >= ControlDiv_y[18] && y < ControlDiv_y[19] )
     			{
/* Move -y */
				up_down -= step_sizey;
     			}
     			if ( y >= ControlDiv_y[19] && y < ControlDiv_y[20] )
     			{
/* Move -z */
				in_out -= step_sizez;
     			}
     			if ( y >= ControlDiv_y[20] && y < ControlDiv_y[21] )
     			{
/* Reset Position */
                		left_right = left_right0;
                		up_down = up_down0;
                		in_out = in_out0;

				ortho_right = ortho_right0;
				ortho_left = ortho_left0;
				ortho_top = ortho_top0;
				ortho_bottom = ortho_bottom0;
     			}
/* These are for the Deformation Keys */
     			if ( y >= ControlDiv_y[24] && y < ControlDiv_y[25] )
     			{
/* Before Turned On */
				After_flag = 0;
				/*amplify_factor = 1.0;
				Amplify_flag = 0;*/
				Before_flag = 1;
				Both_flag = 0;
				/*disp_flag = 0;*/
				strain_flag = 0;
				if ( post_flag )
					stress_flag = 0;
     			}
     			if ( y >= ControlDiv_y[25] && y < ControlDiv_y[26] )
     			{
/* After Turned On */
				After_flag = 1;
				Before_flag = 0;
				Both_flag = 0;
     			}
     			if ( y >= ControlDiv_y[26] && y < ControlDiv_y[27] )
     			{
/* Both Before and After Turned On */
				After_flag = 1;
				/*AppliedForce_flag = 0;*/
				Both_flag = 1;
				Before_flag = 1;
				/*disp_flag = 0;
				AppliedDisp_flag = 0;*/
				/*Material_flag = 0;
				Render_flag = 0;
				Solid_flag = 0;
				strain_flag = 0;
				stress_flag = 0;*/
     			}
     			if ( y >= ControlDiv_y[27] && y < ControlDiv_y[28] )
     			{
/* Amplification increased */
     			    if ( post_flag )
     			    {
			    	After_flag = 1;
				amplify_step = amplify_step0;
				if( amplify_factor < 1.0 - SMALL2 )
                                	amplify_step = .1;
			    	amplify_factor += amplify_step;
			    	Amplify_flag = 1;
			    	/*AppliedForce_flag = 0;
			    	AppliedDisp_flag = 0;*/
				for ( i = 0; i < numnp; ++i )
				{
/* Update Coordinates */
			   	   *(coord + nsd*i) =
					*(coord0+nsd*i) +  (*(U+ndof*i) +
					*(Uz_fib + i)*(*(U+ndof*i+4)))*amplify_factor;
			   	   *(coord + nsd*i+1) =
					*(coord0+nsd*i+1) +  (*(U+ndof*i+1) -
					*(Uz_fib + i)*(*(U+ndof*i+3)))*amplify_factor;
			   	   *(coord + nsd*i+2) =
					*(coord0+nsd*i+2) + *(U+ndof*i+2)*amplify_factor;

			   	   *(coord + nsd*(numnp+i)) =
					*(coord0+nsd*(numnp+i)) +  (*(U+ndof*i) +
					*(Uz_fib + i)*(*(U+ndof*i+4)))*amplify_factor;
			   	   *(coord + nsd*(numnp+i) + 1) =
					*(coord0+nsd*(numnp+i)+1) +  (*(U+ndof*i+1) -
					*(Uz_fib + i)*(*(U+ndof*i+3)))*amplify_factor;
			   	   *(coord + nsd*(numnp+i) + 2) =
					*(coord0+nsd*(numnp+i)+2) + *(U+ndof*i+2)*amplify_factor;
				}
/* Update force graphics vectors */
				check = shnormal_vectors(connecter, coord, norm );
				if(!check) printf( " Problems with shnormal_vectors \n");

                	    	for( i = 0; i < bc.num_force[0]; ++i)
                	    	{
				    fpointx = *(coord+nsd*bc.force[i]);
				    fpointy = *(coord+nsd*bc.force[i] + 1);
				    fpointz = *(coord+nsd*bc.force[i] + 2);
				    force_vec[i].x = fpointx - force_vec0[i].x;
				    force_vec[i].y = fpointy - force_vec0[i].y;
				    force_vec[i].z = fpointz - force_vec0[i].z;
				    force_vec[i].phix = fpointx - force_vec0[i].phix;
				    force_vec[i].phiy = fpointy - force_vec0[i].phiy;
                	    	}
                	    }
     			}

/* These are for the Engineering Analysis Option Keys */
     			if ( y >= ControlDiv_y[30] )
     			{
/* Stresses or displacement Turned On */
			        if( post_flag)
				{
				   After_flag = 1;
				   Before_flag = 0;
				}
				Both_flag = 0;
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 0;
				Render_flag = 0;
				Solid_flag = 1;
/* Stresses */

/* Stresses XX*/
				strain_flag = 0;
				stress_flag = 1;
				color_choice = 10;
     			}
     			if ( y >= ControlDiv_y[31] )
     			{
/* Stresses YY*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 11;
     			}
     			if ( y >= ControlDiv_y[32] )
     			{
/* Stresses XY*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 13;
     			}
     			if ( y >= ControlDiv_y[33] )
     			{
/* Stresses ZX*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 14;
     			}
     			if ( y >= ControlDiv_y[34] )
     			{
/* Stresses YZ*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 15;
     			}
     			if ( y >= ControlDiv_y[35] )
     			{
/* Principal Stresses I*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 16;
     			}
     			if ( y >= ControlDiv_y[38] )
     			{
/* Principal Stresses II*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 17;
     			}
     			if ( y >= ControlDiv_y[39] )
     			{
/* Principal Stresses III*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 18;
     			}

/* Displacement */
     			if ( y >= ControlDiv_y[41] )
     			{
/* Displacement X*/
				angle_flag = 0;
				disp_flag = 1;
				stress_flag = 0;
				color_choice = 19;
     			}
     			if ( y >= ControlDiv_y[43] )
     			{
/* Displacement Y*/
				angle_flag = 0;
				disp_flag = 1;
				stress_flag = 0;
				color_choice = 20;
     			}
     			if ( y >= ControlDiv_y[44] )
     			{
/* Displacement Z*/
				angle_flag = 0;
				disp_flag = 1;
				stress_flag = 0;
				color_choice = 21;
     			}
     		}
		if ( x >= textDiv_xa && x < textDiv_xb )
     		{
     			if ( y >= ControlDiv_y[3] && y <= ControlDiv_y[4] )
     			{
/* Solid Turned Off */
				angle_flag = 0;
				disp_flag = 0;
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 0;
				Render_flag = 0;
				Solid_flag = 0;
				strain_flag = 0;
				stress_flag = 0;
     			}
     			if ( y >= ControlDiv_y[4] && y < ControlDiv_y[5] )
     			{
/* Node ID Turned On increment down */
				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 1;
                		node_choicef += .5;
                		if ( node_choicef > 1.0 )
                		{
                			node_choicef = 0.0;
                		}
                		node_choice -= (int)node_choicef;
                		if ( node_choice < 0 )
                		{
                        		node_choice = 2*numnp-1;
                		}
				color_choice = 31;
				Render_flag = 0;
				strain_flag = 0;
				stress_flag = 0;
     				if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}
     			}
     			if ( y >= ControlDiv_y[5] && y < ControlDiv_y[6] )
     			{
/* Element ID Turned On increment down */
				/*ScreenShot( 350, 700);*/

				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Element_flag = 1;
				Material_flag = 0;
				Node_flag = 0;
                		ele_choicef += .5;
                		if ( ele_choicef > 1.0 )
                		{
                			ele_choicef = 0.0;
                		}
                		ele_choice -= (int)ele_choicef;
                		if ( ele_choice < 0 )
                		{
                        		ele_choice = numel-1;
                		}
				color_choice = 32;
				Render_flag = 0;
				Solid_flag = 1;
				strain_flag = 0;
				stress_flag = 0;
     				if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}
     			}
     			if ( y >= ControlDiv_y[6] && y < ControlDiv_y[7] )
     			{
/* Material Turned On increment down */
				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Element_flag = 0;
				Node_flag = 0;
				Material_flag = 1;
                		matl_choicef += .5;
                		if ( matl_choicef > 1.0 )
                		{
                			matl_choicef = 0.0;
                		}
                		matl_choice -= (int)matl_choicef;
                		if ( matl_choice < 0 )
                		{
                        		matl_choice = nmat-1;
                		}
				color_choice = 30;
				Render_flag = 0;
				Solid_flag = 1;
				strain_flag = 0;
				stress_flag = 0;
     				if ( Both_flag )
     				{
					After_flag = 0;
					Before_flag = 1;
					Both_flag = 0;
     				}
     			}
     			if ( y >= ControlDiv_y[7] && y <= ControlDiv_y[8] )
     			{
/* Fixed Disp Turned Off */
				AppliedDisp_flag = 0;
     			}
     			if ( y >= ControlDiv_y[8] && y < ControlDiv_y[9] )
     			{
/* Applied Force Turned Off */
				AppliedForce_flag = 0;
     			}
     			if ( y >= ControlDiv_y[9] && y < ControlDiv_y[10] )
     			{
/* Axes Turned Off */
				Axes_flag = 0;
     			}
/* These are for the Rotation Keys */
     			if ( y >= ControlDiv_y[12] && y < ControlDiv_y[13] )
     			{
/* Rotate +x */
				xAngle += 5.0;
     			}
     			if ( y >= ControlDiv_y[13] && y < ControlDiv_y[14] )
     			{
/* Rotate +y */
				yAngle += 5.0;
     			}
     			if ( y >= ControlDiv_y[14] && y < ControlDiv_y[15] )
     			{
/* Rotate +z */
				zAngle += 5.0;
     			}
/* These are for the Move Keys */
     			if ( y >= ControlDiv_y[17] && y < ControlDiv_y[18] )
     			{
/* Move +x */
				left_right += step_sizex;
     			}
     			if ( y >= ControlDiv_y[18] && y < ControlDiv_y[19] )
     			{
/* Move +y */
				up_down += step_sizey;
     			}
     			if ( y >= ControlDiv_y[19] && y < ControlDiv_y[20] )
     			{
/* Move +z */
				in_out += step_sizez;
     			}
/* These are for the Deformation Keys */
     			if ( y >= ControlDiv_y[27] && y < ControlDiv_y[28] )
     			{
/* Amplification decreased */
            		    if( post_flag )
            		    {
			    	After_flag = 1;
				amplify_step = amplify_step0;
				if( amplify_factor < 1.0 + amplify_step0 - SMALL2 )
					amplify_step = .1;
			    	amplify_factor -= amplify_step;
			    	/*Amplify_flag = 1;*/
     			    	if ( amplify_factor < 0.0 )
				    amplify_factor = 0.0;
				for ( i = 0; i < numnp; ++i )
				{
/* Update Coordinates */
			   	   *(coord + nsd*i) =
					*(coord0+nsd*i) +  (*(U+ndof*i) +
					*(Uz_fib + i)*(*(U+ndof*i+4)))*amplify_factor;
			   	   *(coord + nsd*i+1) =
					*(coord0+nsd*i+1) +  (*(U+ndof*i+1) -
					*(Uz_fib + i)*(*(U+ndof*i+3)))*amplify_factor;
			   	   *(coord + nsd*i+2) =
					*(coord0+nsd*i+2) + *(U+ndof*i+2)*amplify_factor;

			   	   *(coord + nsd*(numnp+i)) =
					*(coord0+nsd*(numnp+i)) +  (*(U+ndof*i) +
					*(Uz_fib + i)*(*(U+ndof*i+4)))*amplify_factor;
			   	   *(coord + nsd*(numnp+i) + 1) =
					*(coord0+nsd*(numnp+i)+1) +  (*(U+ndof*i+1) -
					*(Uz_fib + i)*(*(U+ndof*i+3)))*amplify_factor;
			   	   *(coord + nsd*(numnp+i) + 2) =
					*(coord0+nsd*(numnp+i)+2) + *(U+ndof*i+2)*amplify_factor;
				}
/* Update force graphics vectors */
				check = shnormal_vectors(connecter, coord, norm );
				if(!check) printf( " Problems with shnormal_vectors \n");

                	    	for( i = 0; i < bc.num_force[0]; ++i)
                	    	{
				    fpointx = *(coord+nsd*bc.force[i]);
				    fpointy = *(coord+nsd*bc.force[i] + 1);
				    fpointz = *(coord+nsd*bc.force[i] + 2);
				    force_vec[i].x = fpointx - force_vec0[i].x;
				    force_vec[i].y = fpointy - force_vec0[i].y;
				    force_vec[i].z = fpointz - force_vec0[i].z;
				    force_vec[i].phix = fpointx - force_vec0[i].phix;
				    force_vec[i].phiy = fpointy - force_vec0[i].phiy;
                	    	}
                	    }
     			}

/* These are for the Engineering Analysis Option Keys */
     			if ( y >= ControlDiv_y[30] )
     			{
/* Strains Turned On */
			        if( post_flag)
				{
				   After_flag = 1;
				   Before_flag = 0;
				}
				Both_flag = 0;
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 0;
				Render_flag = 0;
				Solid_flag = 1;
/* Strains */

/* Strains XX*/
				strain_flag = 1;
				color_choice = 1;
				stress_flag = 0;
     			}
     			if ( y >= ControlDiv_y[31] )
     			{
/* Strains YY*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 2;
     			}
     			if ( y >= ControlDiv_y[32] )
     			{
/* Strains XY*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 4;
     			}
     			if ( y >= ControlDiv_y[33] )
     			{
/* Strains ZX*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 5;
     			}
     			if ( y >= ControlDiv_y[34] )
     			{
/* Strains YZ*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 6;
     			}
     			if ( y >= ControlDiv_y[35] )
     			{
/* Principal Strains I*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 7;
     			}
     			if ( y >= ControlDiv_y[38] )
     			{
/* Principal Strains II*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 8;
     			}
     			if ( y >= ControlDiv_y[39] )
     			{
/* Principal Strains III*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 9;
     			}
/* Angle */
     			if ( y >= ControlDiv_y[41] )
     			{
/* Angle X*/
				angle_flag = 1;
				disp_flag = 0;
				stress_flag = 0;
				color_choice = 22;
     			}
     			if ( y >= ControlDiv_y[43] )
     			{
/* Angle Y*/
				angle_flag = 2;
				disp_flag = 0;
				color_choice = 23;
     			}
     		}
	}

	sprintf( RotateData[0], "%8.2f ", xAngle);
	sprintf( RotateData[1], "%8.2f ", yAngle);
	sprintf( RotateData[2], "%8.2f ", zAngle);

	sprintf( MoveData[0], "%8.2f ", left_right);
	sprintf( MoveData[1], "%8.2f ", up_down);
	sprintf( MoveData[2], "%8.2f ", in_out);

	sprintf( AmplifyData, "%8.2f ", amplify_factor);

    	strcpy(BoxText, "");
	for(i = 29; i < rowdim ; ++i)
	{
		Color_flag[i] = 0;
	}

        sprintf( BoxData[0], " " );
        sprintf( BoxData[1], " " );
        sprintf( BoxData[2], " " );
        sprintf( BoxData[3], " " );
        sprintf( BoxData[4], " " );
        sprintf( BoxData[5], " " );
        sprintf( BoxData[6], " " );
        sprintf( BoxData[7], " " );
        sprintf( BoxData[8], " " );
        sprintf( BoxData[9], " " );
        sprintf( BoxData[10], " " );
        sprintf( BoxData[11], " " );
        sprintf( BoxData[12], " " );
        sprintf( BoxData[13], " " );
        sprintf( BoxData[14], " " );
        sprintf( BoxData[15], " " );
        sprintf( BoxData[16], " " );
        sprintf( BoxData[17], " " );

    	switch (color_choice) {
                case 1:
                        strcpy(BoxText, "strain XX");
			Color_flag[29] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].xx);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].xx);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].xx);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].xx);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].xx);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].xx);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].xx);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].xx);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].xx);
                break;
                case 2:
                        strcpy(BoxText, "strain YY");
			Color_flag[30] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].yy);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].yy);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].yy);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].yy);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].yy);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].yy);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].yy);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].yy);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].yy);
                break;
                case 4:
                        strcpy(BoxText, "strain XY");
			Color_flag[31] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].xy);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].xy);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].xy);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].xy);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].xy);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].xy);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].xy);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].xy);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].xy);
                break;
                case 5:
                        strcpy(BoxText, "strain ZX");
			Color_flag[32] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].zx);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].zx);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].zx);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].zx);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].zx);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].zx);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].zx);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].zx);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].zx);
                break;
                case 6:
                        strcpy(BoxText, "strain YZ");
			Color_flag[33] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].yz);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].yz);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].yz);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].yz);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].yz);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].yz);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].yz);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].yz);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].yz);
                break;
                case 7:
                        strcpy(BoxText, "strain I");
			Color_flag[36] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].I);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].I);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].I);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].I);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].I);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].I);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].I);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].I);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].I);
                break;
                case 8:
                        strcpy(BoxText, "strain II");
			Color_flag[37] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].II);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].II);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].II);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].II);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].II);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].II);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].II);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].II);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].II);
                break;
                case 9:
                        strcpy(BoxText, "strain III");
			Color_flag[38] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", strain_div[8].III);
			sprintf( BoxData[2], "%10.3e ", strain_div[7].III);
			sprintf( BoxData[4], "%10.3e ", strain_div[6].III);
			sprintf( BoxData[6], "%10.3e ", strain_div[5].III);
			sprintf( BoxData[8], "%10.3e ", strain_div[4].III);
			sprintf( BoxData[10], "%10.3e ", strain_div[3].III);
			sprintf( BoxData[12], "%10.3e ", strain_div[2].III);
			sprintf( BoxData[14], "%10.3e ", strain_div[1].III);
			sprintf( BoxData[16], "%10.3e ", strain_div[0].III);
                break;
      		case 10:
    			strcpy(BoxText, "stress XX");
			Color_flag[29] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].xx);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].xx);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].xx);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].xx);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].xx);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].xx);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].xx);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].xx);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].xx);
       		break;
      		case 11:
    			strcpy(BoxText, "stress YY");
			Color_flag[30] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].yy);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].yy);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].yy);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].yy);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].yy);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].yy);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].yy);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].yy);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].yy);
       		break;
      		case 13:
    			strcpy(BoxText, "stress XY");
			Color_flag[31] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].xy);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].xy);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].xy);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].xy);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].xy);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].xy);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].xy);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].xy);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].xy);
       		break;
      		case 14:
    			strcpy(BoxText, "stress ZX");
			Color_flag[32] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].zx);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].zx);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].zx);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].zx);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].zx);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].zx);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].zx);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].zx);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].zx);
       		break;
      		case 15:
    			strcpy(BoxText, "stress YZ");
			Color_flag[33] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].yz);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].yz);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].yz);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].yz);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].yz);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].yz);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].yz);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].yz);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].yz);
       		break;
      		case 16:
    			strcpy(BoxText, "stress I");
			Color_flag[36] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].I);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].I);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].I);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].I);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].I);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].I);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].I);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].I);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].I);
       		break;
      		case 17:
    			strcpy(BoxText, "stress II");
			Color_flag[37] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].II);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].II);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].II);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].II);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].II);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].II);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].II);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].II);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].II);
       		break;
      		case 18:
    			strcpy(BoxText, "stress III");
			Color_flag[38] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", stress_div[8].III);
			sprintf( BoxData[2], "%10.3e ", stress_div[7].III);
			sprintf( BoxData[4], "%10.3e ", stress_div[6].III);
			sprintf( BoxData[6], "%10.3e ", stress_div[5].III);
			sprintf( BoxData[8], "%10.3e ", stress_div[4].III);
			sprintf( BoxData[10], "%10.3e ", stress_div[3].III);
			sprintf( BoxData[12], "%10.3e ", stress_div[2].III);
			sprintf( BoxData[14], "%10.3e ", stress_div[1].III);
			sprintf( BoxData[16], "%10.3e ", stress_div[0].III);
       		break;
                case 19:
                        strcpy(BoxText, "disp X");
			Color_flag[41] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 1;
			sprintf( BoxData[0], "%10.3e ", Ux_div[8]);
			sprintf( BoxData[2], "%10.3e ", Ux_div[7]);
			sprintf( BoxData[4], "%10.3e ", Ux_div[6]);
			sprintf( BoxData[6], "%10.3e ", Ux_div[5]);
			sprintf( BoxData[8], "%10.3e ", Ux_div[4]);
			sprintf( BoxData[10], "%10.3e ", Ux_div[3]);
			sprintf( BoxData[12], "%10.3e ", Ux_div[2]);
			sprintf( BoxData[14], "%10.3e ", Ux_div[1]);
			sprintf( BoxData[16], "%10.3e ", Ux_div[0]);
                break;
                case 20:
                        strcpy(BoxText, "disp Y");
			Color_flag[42] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 2;
			sprintf( BoxData[0], "%10.3e ", Uy_div[8]);
			sprintf( BoxData[2], "%10.3e ", Uy_div[7]);
			sprintf( BoxData[4], "%10.3e ", Uy_div[6]);
			sprintf( BoxData[6], "%10.3e ", Uy_div[5]);
			sprintf( BoxData[8], "%10.3e ", Uy_div[4]);
			sprintf( BoxData[10], "%10.3e ", Uy_div[3]);
			sprintf( BoxData[12], "%10.3e ", Uy_div[2]);
			sprintf( BoxData[14], "%10.3e ", Uy_div[1]);
			sprintf( BoxData[16], "%10.3e ", Uy_div[0]);
                break;
                case 21:
                        strcpy(BoxText, "disp Z");
			Color_flag[43] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 3;
			sprintf( BoxData[0], "%10.3e ", Uz_div[8]);
			sprintf( BoxData[2], "%10.3e ", Uz_div[7]);
			sprintf( BoxData[4], "%10.3e ", Uz_div[6]);
			sprintf( BoxData[6], "%10.3e ", Uz_div[5]);
			sprintf( BoxData[8], "%10.3e ", Uz_div[4]);
			sprintf( BoxData[10], "%10.3e ", Uz_div[3]);
			sprintf( BoxData[12], "%10.3e ", Uz_div[2]);
			sprintf( BoxData[14], "%10.3e ", Uz_div[1]);
			sprintf( BoxData[16], "%10.3e ", Uz_div[0]);
                break;
                case 22:
                        strcpy(BoxText, "Angle X");
			Color_flag[41] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 1;
			disp_flag = 0;
			sprintf( BoxData[0], "%10.3e ", Uphi_x_div[8]);
			sprintf( BoxData[2], "%10.3e ", Uphi_x_div[7]);
			sprintf( BoxData[4], "%10.3e ", Uphi_x_div[6]);
			sprintf( BoxData[6], "%10.3e ", Uphi_x_div[5]);
			sprintf( BoxData[8], "%10.3e ", Uphi_x_div[4]);
			sprintf( BoxData[10], "%10.3e ", Uphi_x_div[3]);
			sprintf( BoxData[12], "%10.3e ", Uphi_x_div[2]);
			sprintf( BoxData[14], "%10.3e ", Uphi_x_div[1]);
			sprintf( BoxData[16], "%10.3e ", Uphi_x_div[0]);
                break;
                case 23:
                        strcpy(BoxText, "Angle Y");
			Color_flag[42] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 2;
			disp_flag = 0;
			sprintf( BoxData[0], "%10.3e ", Uphi_y_div[8]);
			sprintf( BoxData[2], "%10.3e ", Uphi_y_div[7]);
			sprintf( BoxData[4], "%10.3e ", Uphi_y_div[6]);
			sprintf( BoxData[6], "%10.3e ", Uphi_y_div[5]);
			sprintf( BoxData[8], "%10.3e ", Uphi_y_div[4]);
			sprintf( BoxData[10], "%10.3e ", Uphi_y_div[3]);
			sprintf( BoxData[12], "%10.3e ", Uphi_y_div[2]);
			sprintf( BoxData[14], "%10.3e ", Uphi_y_div[1]);
			sprintf( BoxData[16], "%10.3e ", Uphi_y_div[0]);
                break;
                case 30:
    			strcpy(BoxText, "Material");
			Color_flag[7] = 1;
        		input_color_flag = 0;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 0;
			Element_flag = 0;
			Node_flag = 0;
			Material_flag = 1;
			sprintf( BoxData[0], "%4d ", matl_choice);
			strcpy( BoxData[1], "Emod");
			sprintf( BoxData[2], "%10.3e ", matl_crtl[matl_choice].E);
			strcpy( BoxData[3], "Poisson");
			sprintf( BoxData[4], "%10.3e ", matl_crtl[matl_choice].nu);
			strcpy( BoxData[5], "Mass");
			sprintf( BoxData[6], "%10.3e ", matl_crtl[matl_choice].rho);
			strcpy( BoxData[7], "shear fac.");
			sprintf( BoxData[8], "%10.3e ", matl_crtl[matl_choice].shear);
			sprintf( BoxData[14], " " );
			sprintf( BoxData[16], " " );
                break;
                case 31:
    			strcpy(BoxText, "Node");
			Color_flag[3] = 1;
        		input_color_flag = 0;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 0;
			Element_flag = 0;
			Material_flag = 0;
			Node_flag = 1;
			fpointx = *(coord + nsd*node_choice);
			fpointy = *(coord + nsd*node_choice + 1);
			fpointz = *(coord + nsd*node_choice + 2);
			if(!After_flag)
			{
			        fpointx = *(coord0 + nsd*node_choice);
			        fpointy = *(coord0 + nsd*node_choice + 1);
			        fpointz = *(coord0 + nsd*node_choice + 2);
			}
			sprintf( BoxData[0], "%4d ", node_choice);
			strcpy( BoxData[2], "coord x");
			sprintf( BoxData[4], "%10.3e ", fpointx);
			strcpy( BoxData[6], "coord y");
			sprintf( BoxData[8], "%10.3e ", fpointy);
			strcpy( BoxData[10], "coord z");
			sprintf( BoxData[12], "%10.3e ", fpointz);
			sprintf( BoxData[14], " " );
			sprintf( BoxData[16], " " );
                break;
                case 32:
    			strcpy(BoxText, "Element");
			Color_flag[4] = 1;
        		input_color_flag = 0;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 0;
			Material_flag = 0;
			Node_flag = 0;
			Element_flag = 1;
			sprintf( BoxData[0], "%4d ", ele_choice);
    			strcpy( BoxData[2], "Material");
			sprintf( BoxData[4], "%4d ", *(el_matl_color+ele_choice));
    			strcpy( BoxData[6], "Connect");
			dum1 = *(connecter + npell*ele_choice);
			dum2 = *(connecter + npell*ele_choice+1);
			sprintf( BoxData[8], "%4d,%4d ",dum1, dum2);
			dum1 = *(connecter + npell*ele_choice+2);
			dum2 = *(connecter + npell*ele_choice+3);
			sprintf( BoxData[10], "%4d,%4d ",dum1, dum2);
			sprintf( BoxData[12], " " );
			sprintf( BoxData[14], " " );
			sprintf( BoxData[16], " " );
                break;
        }

        input_color_flag = 0;

/* If there is a post file, then turn the input_color_flag on so that the before
   mesh will be drawn in pink.  If there is no post file, turn on the
   input_color_flag for every case except when stress analysis or material, element
   or node is selected.
 */

        if( color_choice < 10)
             input_color_flag = 1;
        if( color_choice > 15 && color_choice < 19)
             input_color_flag = 1;
        if( post_flag > 0 && color_choice < 30)
             input_color_flag = 1;

	Color_flag[2] = Solid_flag;
/*
	Color_flag[3] = Perspective_flag;
	Color_flag[4] = Render_flag;
*/
	Color_flag[3] = Node_flag;
	Color_flag[4] = Element_flag;
	Color_flag[5] = Material_flag;
	Color_flag[6] = AppliedDisp_flag;
	Color_flag[7] = AppliedForce_flag;
	Color_flag[8] = Axes_flag;

	Color_flag[23] = Before_flag;
	Color_flag[24] = After_flag;
	Color_flag[25] = Both_flag;
	Color_flag[26] = Amplify_flag;

        if(!post_flag) After_flag = 0;
        if(!input_flag) Before_flag = 0;

	glutPostWindowRedisplay(ControlWindow);
	glutPostWindowRedisplay(MeshWindow);
}


