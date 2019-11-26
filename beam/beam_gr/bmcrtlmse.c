/*
    This program contains the control mouse routine for the FEM GUI
    for beam elements.
  
                        Last Update 2/22/02

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

/****** FEA globals ******/

extern int dof, nmat, numnp, numel;
extern double *coord, *coord0;
extern double *U;
extern int *connecter;
extern BOUND bc;
extern int *el_type;
extern XYZPhiF *force_vec, *force_vec0;
extern QYQZ *dist_load_vec0;
extern XYZF_GR *dist_load_vec;

/* Global variables for the mesh color and nodal data */

extern int *el_matl_color;
extern MATL *matl_crtl;

/****** EXTERNAL VARIABLES ********/

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
extern int Dist_Load_flag, Perspective_flag, Render_flag, AppliedDisp_flag,
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
extern double Ux_div[boxnumber+1], Uy_div[boxnumber+1], Uz_div[boxnumber+1];
extern double Uphi_x_div[boxnumber+1], Uphi_y_div[boxnumber+1], Uphi_z_div[boxnumber+1];
extern MDIM moment_div[boxnumber+1], curve_div[boxnumber+1];
extern SDIM stress_div[boxnumber+1], strain_div[boxnumber+1];

void ScreenShot( int , int );

int rotater( double *, double *, double *);

void bmControlMouse(int button, int state, int x, int y)
{
	int i, k, check, node0, node1, dum1, dum2;
	double fpointx, fpointy, fpointz;
	double vec_in[3], vec_out[3], coord_el[3*npel];
  	if (button == GLUT_LEFT_BUTTON)
  	{
		if ( x < textDiv_xa )
     		{

/* These are for the View Option Keys */

     			if ( y >= ControlDiv_y[3] && y < ControlDiv_y[4] )
     			{
/* Distributed Load Turned On */
				Dist_Load_flag = 1;
     			}
     			if ( y >= ControlDiv_y[4] && y < ControlDiv_y[5] )
     			{
/* Node ID Turned On increment up */
				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Dist_Load_flag = 0;
				Element_flag = 0;
				Material_flag = 0;
				Node_flag = 1;
                		node_choicef += .5;
                		if ( node_choicef > 1.0 )
                		{
                			node_choicef = 0.0;
                		}
                		node_choice += (int)node_choicef;
                		if ( node_choice > numnp - 1 )
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
				Dist_Load_flag = 0;
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
				Dist_Load_flag = 0;
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
			   	    *(coord + nsd*i) = *(coord0+nsd*i) +
					*(U+ndof*i)*amplify_factor;
			   	    *(coord + nsd*i+1) = *(coord0+nsd*i+1) +
					*(U+ndof*i+1)*amplify_factor;
			   	    *(coord + nsd*i+2) = *(coord0+nsd*i+2) +
					*(U+ndof*i+2)*amplify_factor;
                        	}

/* Update force graphics vectors */	
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
			   	    force_vec[i].phiz = fpointz - force_vec0[i].phiz;
				}
/* Update distributed load graphics vectors */	
            			for( k = 0; k < bc.num_dist_load[0]; ++k)
            			{
                			node0 = *(connecter+bc.dist_load[k]*npel);
                			node1 = *(connecter+bc.dist_load[k]*npel+1);

                			*(coord_el)=*(coord+nsd*node0);
                			*(coord_el+1)=*(coord+nsd*node0+1);
                			*(coord_el+2)=*(coord+nsd*node0+2);

                			*(coord_el+3)=*(coord+nsd*node1);
                			*(coord_el+4)=*(coord+nsd*node1+1);
                			*(coord_el+5)=*(coord+nsd*node1+2);

			        	*(vec_in) =  0.0;
                			*(vec_in+1) =  dist_load_vec0[k].qy;
                			*(vec_in+2) =  dist_load_vec0[k].qz;

                			check = rotater(coord_el, vec_in, vec_out);
                			if(!check) printf( " Problems with rotater \n");

                			dist_load_vec[k].x = *(vec_out);
                			dist_load_vec[k].y = *(vec_out+1);
                			dist_load_vec[k].z = *(vec_out+2);
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
/* Moments and Stresses */

/* Moments XX*/
				strain_flag = 0;
				stress_flag = 1;
				color_choice = 10;
     			}
     			if ( y >= ControlDiv_y[31] )
     			{
/* Moments YY*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 11;
     			}
     			if ( y >= ControlDiv_y[32] )
     			{
/* Moments ZZ*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 12;
     			}
     			if ( y >= ControlDiv_y[33] )
     			{
/* Stresses XX*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 13;
     			}

/* Displacement */
     			if ( y >= ControlDiv_y[36] )
     			{
/* Displacement X*/
				angle_flag = 0;
				disp_flag = 1;
				stress_flag = 0;
				color_choice = 19;
     			}
     			if ( y >= ControlDiv_y[38] )
     			{
/* Displacement Y*/
				angle_flag = 0;
				disp_flag = 1;
				stress_flag = 0;
				color_choice = 20;
     			}
     			if ( y >= ControlDiv_y[39] )
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
/* Distributed Load Turned Off */
				Dist_Load_flag = 0;
     			}
     			if ( y >= ControlDiv_y[4] && y < ControlDiv_y[5] )
     			{
/* Node ID Turned On increment down */
				AppliedForce_flag = 0;
				angle_flag = 0;
				disp_flag = 0;
				AppliedDisp_flag = 0;
				Dist_Load_flag = 0;
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
                        		node_choice = numnp-1;
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
				Dist_Load_flag = 0;
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
				Dist_Load_flag = 0;
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
			   	    *(coord + nsd*i) = *(coord0+nsd*i) +
					*(U+ndof*i)*amplify_factor;
			   	    *(coord + nsd*i+1) = *(coord0+nsd*i+1) +
					*(U+ndof*i+1)*amplify_factor;
			   	    *(coord + nsd*i+2) = *(coord0+nsd*i+2) +
					*(U+ndof*i+2)*amplify_factor;
                        	}

/* Update force graphics vectors */	
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
			   	    force_vec[i].phiz = fpointz - force_vec0[i].phiz;
				}
/* Update distributed load graphics vectors */	
            			for( k = 0; k < bc.num_dist_load[0]; ++k)
            			{
                			node0 = *(connecter+bc.dist_load[k]*npel);
                			node1 = *(connecter+bc.dist_load[k]*npel+1);

                			*(coord_el)=*(coord+nsd*node0);
                			*(coord_el+1)=*(coord+nsd*node0+1);
                			*(coord_el+2)=*(coord+nsd*node0+2);

                			*(coord_el+3)=*(coord+nsd*node1);
                			*(coord_el+4)=*(coord+nsd*node1+1);
                			*(coord_el+5)=*(coord+nsd*node1+2);

			        	*(vec_in) =  0.0;
                			*(vec_in+1) =  dist_load_vec0[k].qy;
                			*(vec_in+2) =  dist_load_vec0[k].qz;

                			check = rotater(coord_el, vec_in, vec_out);
                			if(!check) printf( " Problems with rotater \n");

                			dist_load_vec[k].x = *(vec_out);
                			dist_load_vec[k].y = *(vec_out+1);
                			dist_load_vec[k].z = *(vec_out+2);
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
/* Curvatures and Strains */

/* Curvature XX*/
				strain_flag = 1;
				color_choice = 1;
				stress_flag = 0;
     			}
     			if ( y >= ControlDiv_y[31] )
     			{
/* Curvature YY*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 2;
     			}
     			if ( y >= ControlDiv_y[32] )
     			{
/* Curvature ZZ*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 3;
     			}
     			if ( y >= ControlDiv_y[33] )
     			{
/* Strains XX*/
				angle_flag = 0;
				disp_flag = 0;
				color_choice = 4;
     			}
/* Angle */
     			if ( y >= ControlDiv_y[36] )
     			{
/* Angle X*/
				angle_flag = 1;
				disp_flag = 0;
				stress_flag = 0;
				color_choice = 22;
     			}
     			if ( y >= ControlDiv_y[38] )
     			{
/* Angle Y*/
				angle_flag = 2;
				disp_flag = 0;
				color_choice = 23;
     			}
     			if ( y >= ControlDiv_y[39] )
     			{
/* Angle Z*/
				angle_flag = 3;
				disp_flag = 0;
				color_choice = 24;
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
                        strcpy(BoxText, "curve XX");
			Color_flag[29] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", curve_div[8].xx);
			sprintf( BoxData[2], "%10.3e ", curve_div[7].xx);
			sprintf( BoxData[4], "%10.3e ", curve_div[6].xx);
			sprintf( BoxData[6], "%10.3e ", curve_div[5].xx);
			sprintf( BoxData[8], "%10.3e ", curve_div[4].xx);
			sprintf( BoxData[10], "%10.3e ", curve_div[3].xx);
			sprintf( BoxData[12], "%10.3e ", curve_div[2].xx);
			sprintf( BoxData[14], "%10.3e ", curve_div[1].xx);
			sprintf( BoxData[16], "%10.3e ", curve_div[0].xx);
                break;
                case 2:
                        strcpy(BoxText, "curve YY");
			Color_flag[30] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", curve_div[8].yy);
			sprintf( BoxData[2], "%10.3e ", curve_div[7].yy);
			sprintf( BoxData[4], "%10.3e ", curve_div[6].yy);
			sprintf( BoxData[6], "%10.3e ", curve_div[5].yy);
			sprintf( BoxData[8], "%10.3e ", curve_div[4].yy);
			sprintf( BoxData[10], "%10.3e ", curve_div[3].yy);
			sprintf( BoxData[12], "%10.3e ", curve_div[2].yy);
			sprintf( BoxData[14], "%10.3e ", curve_div[1].yy);
			sprintf( BoxData[16], "%10.3e ", curve_div[0].yy);
                break;
                case 3:
                        strcpy(BoxText, "curve ZZ");
			Color_flag[31] = 1;
			strain_flag = 1;
			stress_flag = 0;
			sprintf( BoxData[0], "%10.3e ", curve_div[8].zz);
			sprintf( BoxData[2], "%10.3e ", curve_div[7].zz);
			sprintf( BoxData[4], "%10.3e ", curve_div[6].zz);
			sprintf( BoxData[6], "%10.3e ", curve_div[5].zz);
			sprintf( BoxData[8], "%10.3e ", curve_div[4].zz);
			sprintf( BoxData[10], "%10.3e ", curve_div[3].zz);
			sprintf( BoxData[12], "%10.3e ", curve_div[2].zz);
			sprintf( BoxData[14], "%10.3e ", curve_div[1].zz);
			sprintf( BoxData[16], "%10.3e ", curve_div[0].zz);
                break;
                case 4:
                        strcpy(BoxText, "strain XX");
			Color_flag[33] = 1;
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
      		case 10:
    			strcpy(BoxText, "moment XX");
			Color_flag[29] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", moment_div[8].xx);
			sprintf( BoxData[2], "%10.3e ", moment_div[7].xx);
			sprintf( BoxData[4], "%10.3e ", moment_div[6].xx);
			sprintf( BoxData[6], "%10.3e ", moment_div[5].xx);
			sprintf( BoxData[8], "%10.3e ", moment_div[4].xx);
			sprintf( BoxData[10], "%10.3e ", moment_div[3].xx);
			sprintf( BoxData[12], "%10.3e ", moment_div[2].xx);
			sprintf( BoxData[14], "%10.3e ", moment_div[1].xx);
			sprintf( BoxData[16], "%10.3e ", moment_div[0].xx);
       		break;
      		case 11:
    			strcpy(BoxText, "moment YY");
			Color_flag[30] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", moment_div[8].yy);
			sprintf( BoxData[2], "%10.3e ", moment_div[7].yy);
			sprintf( BoxData[4], "%10.3e ", moment_div[6].yy);
			sprintf( BoxData[6], "%10.3e ", moment_div[5].yy);
			sprintf( BoxData[8], "%10.3e ", moment_div[4].yy);
			sprintf( BoxData[10], "%10.3e ", moment_div[3].yy);
			sprintf( BoxData[12], "%10.3e ", moment_div[2].yy);
			sprintf( BoxData[14], "%10.3e ", moment_div[1].yy);
			sprintf( BoxData[16], "%10.3e ", moment_div[0].yy);
       		break;
      		case 12:
    			strcpy(BoxText, "moment ZZ");
			Color_flag[31] = 1;
			strain_flag = 0;
			stress_flag = 1;
			sprintf( BoxData[0], "%10.3e ", moment_div[8].zz);
			sprintf( BoxData[2], "%10.3e ", moment_div[7].zz);
			sprintf( BoxData[4], "%10.3e ", moment_div[6].zz);
			sprintf( BoxData[6], "%10.3e ", moment_div[5].zz);
			sprintf( BoxData[8], "%10.3e ", moment_div[4].zz);
			sprintf( BoxData[10], "%10.3e ", moment_div[3].zz);
			sprintf( BoxData[12], "%10.3e ", moment_div[2].zz);
			sprintf( BoxData[14], "%10.3e ", moment_div[1].zz);
			sprintf( BoxData[16], "%10.3e ", moment_div[0].zz);
       		break;
      		case 13:
    			strcpy(BoxText, "stress XX");
			Color_flag[33] = 1;
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
                case 19:
                        strcpy(BoxText, "disp X");
			Color_flag[36] = 1;
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
			Color_flag[37] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 1;
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
			Color_flag[38] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 0;
			disp_flag = 1;
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
			Color_flag[36] = 1;
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
			Color_flag[37] = 1;
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
                case 24:
                        strcpy(BoxText, "Angle Z");
			Color_flag[38] = 1;
			strain_flag = 0;
			stress_flag = 0;
			angle_flag = 3;
			disp_flag = 0;
			sprintf( BoxData[0], "%10.3e ", Uphi_z_div[8]);
			sprintf( BoxData[2], "%10.3e ", Uphi_z_div[7]);
			sprintf( BoxData[4], "%10.3e ", Uphi_z_div[6]);
			sprintf( BoxData[6], "%10.3e ", Uphi_z_div[5]);
			sprintf( BoxData[8], "%10.3e ", Uphi_z_div[4]);
			sprintf( BoxData[10], "%10.3e ", Uphi_z_div[3]);
			sprintf( BoxData[12], "%10.3e ", Uphi_z_div[2]);
			sprintf( BoxData[14], "%10.3e ", Uphi_z_div[1]);
			sprintf( BoxData[16], "%10.3e ", Uphi_z_div[0]);
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
    			strcpy( BoxData[7], "Area");
			sprintf( BoxData[8], "%10.3e ", matl_crtl[matl_choice].area);
    			strcpy( BoxData[9], "Iy");
			sprintf( BoxData[10], "%10.3e ", matl_crtl[matl_choice].Iy);
    			strcpy( BoxData[11], "Iz");
			sprintf( BoxData[12], "%10.3e ", matl_crtl[matl_choice].Iz);
			sprintf( BoxData[13], " ");
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
    			strcpy( BoxData[6], "type");
			if( *(el_type+ele_choice) > 1)
			{
    				strcpy( BoxData[8], "  beam");
			}
			else
			{
    				strcpy( BoxData[8], "  truss");
			}
    			strcpy( BoxData[10], "Connect");
			dum1 = *(connecter + npel*ele_choice);
			dum2 = *(connecter + npel*ele_choice+1);
			sprintf( BoxData[12], "%4d,%4d ",dum1, dum2);
			sprintf( BoxData[13], " " );
			sprintf( BoxData[14], " " );
			sprintf( BoxData[15], " " );
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

	Color_flag[2] = Dist_Load_flag;
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


