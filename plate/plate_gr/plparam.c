/*
    This program calculates and writes the parameters for
    the FEM GUI for plate elements.
  
   			Last Update 6/9/00

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
 */
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../plate/plconst.h"
#include "../plate/plstruct.h"
#include "plstrcgr.h"
#include "../../common_gr/control.h"

#define init_far0      -2.0

extern int nmat, numnp, numel, dof;
extern double step_sizex, step_sizey, step_sizez;
extern double left, right, top, bottom, near, far, fscale;
extern int control_height, control_width, mesh_height, mesh_width;
extern double ortho_left, ortho_right, ortho_top, ortho_bottom,
	ortho_left0, ortho_right0, ortho_top0, ortho_bottom0;
extern double left_right, up_down, in_out, left_right0, up_down0, in_out0;
extern double AxisMax_x, AxisMax_y, AxisMax_z,
	AxisMin_x, AxisMin_y, AxisMin_z,
	IAxisMin_x, IAxisMin_y, IAxisMin_z;
extern double AxisLength_x, AxisLength_y, AxisLength_z,
	AxisLength_max;
extern double AxisPoint_step;
extern double amplify_step0;

extern double init_right, init_left, init_top,
	init_bottom, init_near, init_far, true_far, dim_max;
extern MDIM del_moment, del_curve, max_moment, min_moment,
	max_curve, min_curve;
extern SDIM del_stress, del_strain, max_stress, min_stress,
	max_strain, min_strain;
extern double max_Uphi_x, min_Uphi_x, del_Uphi_x, max_Uphi_y, min_Uphi_y, del_Uphi_y,
       	max_Uz, min_Uz, del_Uz, absolute_max_U;

int plparameter( double *coord, MDIM *curve_node, MDIM *moment_node, SDIM *strain_node,
	SDIM *stress_node, double *U )
{
        int i, j, check;
        int node_Uz_max, node_Uz_min, node_Uphi_x_max, node_Uphi_x_min,
		node_Uphi_y_max, node_Uphi_y_min;
        IMDIM max_moment_node, min_moment_node, max_curve_node, min_curve_node;
        ISDIM max_stress_node, min_stress_node, max_strain_node, min_strain_node;
	FILE *pldata;

/*   pldata contains all the parameters and extreme values
*/
	pldata = fopen( "plview.dat","w" );

/* Initialize parameters */

	step_sizex = .1; step_sizey = .1; step_sizez = .1;
        init_right = - BIG; init_left = BIG;
	init_top = - BIG; init_bottom = BIG;
	init_near = 0.0; init_far = 0.0; true_far = 0.0;
	max_Uz = -BIG; min_Uz = BIG;
	max_Uphi_x = - BIG; min_Uphi_x = BIG;
	max_Uphi_y = - BIG; min_Uphi_y = BIG;

        node_Uz_max = 0; node_Uz_min = IBIG;
        node_Uphi_x_max = 0; node_Uphi_x_min = IBIG;
        node_Uphi_y_max = 0; node_Uphi_y_min = IBIG;

        max_curve_node.xx = 0; min_curve_node.xx = IBIG;
        max_curve_node.yy = 0; min_curve_node.yy = IBIG;
        max_curve_node.xy = 0; min_curve_node.xy = IBIG;
        max_strain_node.zx = 0; min_strain_node.zx = IBIG;
        max_strain_node.yz = 0; min_strain_node.yz = IBIG;
        max_curve_node.I = 0; min_curve_node.I = IBIG;
        max_curve_node.II = 0; min_curve_node.II = IBIG;

        max_moment_node.xx = 0; min_moment_node.xx = IBIG;
        max_moment_node.yy = 0; min_moment_node.yy = IBIG;
        max_moment_node.xy = 0; min_moment_node.xy = IBIG;
        max_stress_node.zx = 0; min_stress_node.zx = IBIG;
        max_stress_node.yz = 0; min_stress_node.yz = IBIG;
        max_moment_node.I = 0; min_moment_node.I = IBIG;
        max_moment_node.II = 0; min_moment_node.II = IBIG;

/* Initialize for largest and smallest curvatures and strains */

        max_curve.xx = - BIG; min_curve.xx = BIG;
        max_curve.yy = - BIG; min_curve.yy = BIG;
        max_curve.xy = - BIG; min_curve.xy = BIG;
        max_strain.zx = - BIG; min_strain.zx = BIG;
        max_strain.yz = - BIG; min_strain.yz = BIG;
        max_curve.I = - BIG; min_curve.I = BIG;
        max_curve.II = - BIG; min_curve.II = BIG;

/* Initialize largest and smallest moments and stresses */

        max_moment.xx = - BIG; min_moment.xx = BIG;
        max_moment.yy = - BIG; min_moment.yy = BIG;
        max_moment.xy = - BIG; min_moment.xy = BIG;
        max_stress.zx = - BIG; min_stress.zx = BIG;
        max_stress.yz = - BIG; min_stress.yz = BIG;
        max_moment.I = - BIG; min_moment.I = BIG;
        max_moment.II = - BIG; min_moment.II = BIG;

/* Search for extreme values */
 
/* Search for extreme values of displacement, angle, and nodal point */

        for( i = 0; i < numnp; ++i )
        {

/* Search for extreme nodal coordinates for parameters when object is
   viewed orthographically */

                if( init_right < *(coord+nsd*i))
                        init_right=*(coord+nsd*i);

                if( init_left > *(coord+nsd*i))
                        init_left=*(coord+nsd*i);

                if( init_top < *(coord+nsd*i+1))
                        init_top=*(coord+nsd*i+1);

                if( init_bottom > *(coord+nsd*i+1))
                        init_bottom=*(coord+nsd*i+1);

/* Search for extreme nodal displacements and angles */

                if( max_Uz < *(U+ndof*i))
		{
                        max_Uz=*(U+ndof*i);
			node_Uz_max = i;
		}
                if( min_Uz > *(U+ndof*i))
		{
                        min_Uz=*(U+ndof*i);
			node_Uz_min = i;
		}
                if( max_Uphi_x < *(U+ndof*i+1))
		{
                        max_Uphi_x=*(U+ndof*i+1);
			node_Uphi_x_max = i;
		}
                if( min_Uphi_x > *(U+ndof*i+1))
		{
                        min_Uphi_x=*(U+ndof*i+1);
			node_Uphi_x_min = i;
		}
                if( max_Uphi_y < *(U+ndof*i+2))
		{
                        max_Uphi_y=*(U+ndof*i+2);
			node_Uphi_y_max = i;
		}
                if( min_Uphi_y > *(U+ndof*i+2))
		{
                        min_Uphi_y=*(U+ndof*i+2);
			node_Uphi_y_min = i;
		}
        }

/* Search for largest absolute value of displacement U */

	absolute_max_U = fabs(min_Uz);
        if( absolute_max_U < fabs(max_Uz))
		absolute_max_U = fabs(max_Uz);

        if( init_far > true_far)
		init_far=true_far;


        for( i = 0; i < numnp; ++i )
        {

/* Find extreme curvatures and strains */

                if( max_curve.xx < curve_node[i].xx )
		{
                        max_curve.xx = curve_node[i].xx;
			max_curve_node.xx = i;
		}
                if( min_curve.xx > curve_node[i].xx )
		{
                        min_curve.xx = curve_node[i].xx;
			min_curve_node.xx = i;
		}
                if( max_curve.yy < curve_node[i].yy )
		{
                        max_curve.yy = curve_node[i].yy;
			max_curve_node.yy = i;
		}
                if( min_curve.yy > curve_node[i].yy )
		{
                        min_curve.yy = curve_node[i].yy;
			min_curve_node.yy = i;
		}
                if( max_curve.xy < curve_node[i].xy )
		{
                        max_curve.xy = curve_node[i].xy;
			max_curve_node.xy = i;
		}
                if( min_curve.xy > curve_node[i].xy )
		{
                        min_curve.xy = curve_node[i].xy;
			min_curve_node.xy = i;
		}
                if( max_strain.zx < strain_node[i].zx )
		{
                        max_strain.zx = strain_node[i].zx;
			max_strain_node.zx = i;
		}
                if( min_strain.zx > strain_node[i].zx )
		{
                        min_strain.zx = strain_node[i].zx;
			min_strain_node.zx = i;
		}
                if( max_strain.yz < strain_node[i].yz )
		{
                        max_strain.yz = strain_node[i].yz;
			max_strain_node.yz = i;
		}
                if( min_strain.yz > strain_node[i].yz )
		{
                        min_strain.yz = strain_node[i].yz;
			min_strain_node.yz = i;
		}
                if( max_curve.I < curve_node[i].I )
		{
                        max_curve.I = curve_node[i].I;
			max_curve_node.I = i;
		}
                if( min_curve.I > curve_node[i].I )
		{
                        min_curve.I = curve_node[i].I;
			min_curve_node.I = i;
		}
                if( max_curve.II < curve_node[i].II )
		{
                        max_curve.II = curve_node[i].II;
			max_curve_node.II = i;
		}
                if( min_curve.II > curve_node[i].II )
		{
                        min_curve.II = curve_node[i].II;
			min_curve_node.II = i;
		}
/* Find extreme moments and stresses */

                if( max_moment.xx < moment_node[i].xx )
		{
                        max_moment.xx = moment_node[i].xx;
			max_moment_node.xx = i;
		}
                if( min_moment.xx > moment_node[i].xx )
		{
                        min_moment.xx = moment_node[i].xx;
			min_moment_node.xx = i;
		}
                if( max_moment.yy < moment_node[i].yy )
		{
                        max_moment.yy = moment_node[i].yy;
			max_moment_node.yy = i;
		}
                if( min_moment.yy > moment_node[i].yy )
		{
                        min_moment.yy = moment_node[i].yy;
			min_moment_node.yy = i;
		}
                if( max_moment.xy < moment_node[i].xy )
		{
                        max_moment.xy = moment_node[i].xy;
			max_moment_node.xy = i;
		}
                if( min_moment.xy > moment_node[i].xy )
		{
                        min_moment.xy = moment_node[i].xy;
			min_moment_node.xy = i;
		}
                if( max_stress.zx < stress_node[i].zx )
		{
                        max_stress.zx = stress_node[i].zx;
			max_stress_node.zx = i;
		}
                if( min_stress.zx > stress_node[i].zx )
		{
                        min_stress.zx = stress_node[i].zx;
			min_stress_node.zx = i;
		}
                if( max_stress.yz < stress_node[i].yz )
		{
                        max_stress.yz = stress_node[i].yz;
			max_stress_node.yz = i;
		}
                if( min_stress.yz > stress_node[i].yz )
		{
                        min_stress.yz = stress_node[i].yz;
			min_stress_node.yz = i;
		}
                if( max_moment.I < moment_node[i].I )
		{
                        max_moment.I = moment_node[i].I;
			max_moment_node.I = i;
		}
                if( min_moment.I > moment_node[i].I )
		{
                        min_moment.I = moment_node[i].I;
			min_moment_node.I = i;
		}
                if( max_moment.II < moment_node[i].II )
		{
                        max_moment.II = moment_node[i].II;
			max_moment_node.II = i;
		}
                if( min_moment.II > moment_node[i].II )
		{
                        min_moment.II = moment_node[i].II;
			min_moment_node.II = i;
		}
        }

/* Set the axes parameters */

	AxisMax_x = 1.2*init_right + 1.0;
	AxisMax_y = 1.2*init_top + 1.0;
	AxisMin_x = 1.2*init_left - 1.0;
	AxisMin_y = 1.2*init_bottom - 1.0;

	if( AxisMax_x < 0.0 )
		AxisMax_x = -.1*AxisMin_x;
	if( AxisMax_y < 0.0 )
		AxisMax_y = -.1*AxisMin_y;
	if( AxisMin_x > 0.0 )
		AxisMin_x = -.1*AxisMax_x;
	if( AxisMin_y > 0.0 )
		AxisMin_y = -.1*AxisMax_y;

	IAxisMin_x = (double)(int)AxisMin_x; 
	IAxisMin_y = (double)(int)AxisMin_y;

	AxisLength_x = (int)(AxisMax_x - AxisMin_x);
	AxisLength_y = (int)(AxisMax_y - AxisMin_y);
	AxisLength_max = AxisLength_x;
	if( AxisLength_max < AxisLength_y )
		AxisLength_max = AxisLength_y;

	AxisPoint_step = .1;
	if( AxisLength_max > 1.0 )
		AxisPoint_step = 1.0;
	if( AxisLength_max > 10.0 )
		AxisPoint_step = 10.0;
	if( AxisLength_max > 100.0 )
		AxisPoint_step = 100.0;
	if( AxisLength_max > 1000.0 )
		AxisPoint_step = 1000.0;
	if( AxisLength_max > 10000.0 )
		AxisPoint_step = 10000.0;

	AxisLength_max *= .05;

	AxisMax_z = AxisLength_max;
	AxisMin_z = -AxisLength_max;

	IAxisMin_z = (double)(int)AxisMin_z;

/* Determine amplification step size */

	amplify_step0 = .5*AxisLength_max/(absolute_max_U+SMALL);
	/*printf(" amplify_step0 %10.5e\n",amplify_step0);*/

/* Calculate orthographic viewport parameters */

        right = init_right + (init_right - init_left) / 10.0;
        left = init_left - (init_right - init_left) / 10.0;
        top = init_top + (init_top - init_bottom) / 10.0;
        bottom = init_bottom - (init_top - init_bottom) / 10.0;
	near = init_near + (init_near - init_far) / 10.0;
	far = init_far - (init_near - init_far) / 10.0;

	dim_max = right - left;
	if( dim_max < top - bottom )
		dim_max = top - bottom;

	ortho_right0 = left + dim_max + 1.0;
	ortho_left0 = left - 1.0;
	ortho_top0 = bottom + dim_max + 1.0;
	ortho_bottom0 = bottom - 1.0;

	ortho_right = ortho_right0;
	ortho_left = ortho_left0;
	ortho_top = ortho_top0;
	ortho_bottom = ortho_bottom0;

/* Set the Viewer parameters */

	step_sizex = (right - left) / 20.0;
	step_sizey = (top - bottom) / 20.0;
	step_sizez = ( step_sizex + step_sizey ) / 5.0;

	left_right0 = -(left + right) / 2.0;
	up_down0 = - (top + bottom ) / 2.0;
	/*in_out  = (far + near ) / 2.0 - 5.0;*/
	in_out0 = far - 20.0*AxisLength_max;
	left_right = left_right0;
	up_down = up_down0;
	in_out = in_out0;

    	mesh_width=mesh_width0;
    	mesh_height=mesh_height0;

/* Print the above data in the file "plview.dat" */

	fprintf( pldata, "                            node\n");
	fprintf( pldata, "                          min  max       min            max\n");
	fprintf( pldata,"displacement Uz        %5d %5d   %14.6e %14.6e\n", node_Uz_min,
		node_Uz_max, min_Uz, max_Uz);
	fprintf( pldata,"angle phi x            %5d %5d   %14.6e %14.6e\n", node_Uphi_x_min,
		node_Uphi_x_max, min_Uphi_x, max_Uphi_x);
	fprintf( pldata,"angle phi y            %5d %5d   %14.6e %14.6e\n", node_Uphi_y_min,
		node_Uphi_y_max, min_Uphi_y, max_Uphi_y);
	fprintf( pldata,"\n");
	fprintf( pldata, "                            node\n");
	fprintf( pldata, "                        min       max         min           max\n");
	fprintf( pldata,"moment xx            %5d     %5d   %14.6e %14.6e\n", min_moment_node.xx,
		max_moment_node.xx, min_moment.xx, max_moment.xx);
	fprintf( pldata,"moment yy            %5d     %5d   %14.6e %14.6e\n", min_moment_node.yy,
		max_moment_node.yy, min_moment.yy, max_moment.yy);
	fprintf( pldata,"moment xy            %5d     %5d   %14.6e %14.6e\n", min_moment_node.xy,
		max_moment_node.xy, min_moment.xy, max_moment.xy);
	fprintf( pldata,"stress zx            %5d     %5d   %14.6e %14.6e\n", min_stress_node.zx,
		max_stress_node.zx, min_stress.zx, max_stress.zx);
	fprintf( pldata,"stress yz            %5d     %5d   %14.6e %14.6e\n", min_stress_node.yz,
		max_stress_node.yz, min_stress.yz, max_stress.yz);
	fprintf( pldata,"moment I             %5d     %5d   %14.6e %14.6e\n", min_moment_node.I,
		max_moment_node.I, min_moment.I, max_moment.I);
	fprintf( pldata,"moment II            %5d     %5d   %14.6e %14.6e\n", min_moment_node.II,
		max_moment_node.II, min_moment.II, max_moment.II);
	fprintf( pldata,"\n");
	fprintf( pldata,"curve xx             %5d     %5d   %14.6e %14.6e\n", min_curve_node.xx,
		max_curve_node.xx, min_curve.xx, max_curve.xx);
	fprintf( pldata,"curve yy             %5d     %5d   %14.6e %14.6e\n", min_curve_node.yy,
		max_curve_node.yy, min_curve.yy, max_curve.yy);
	fprintf( pldata,"curve xy             %5d     %5d   %14.6e %14.6e\n", min_curve_node.xy,
		max_curve_node.xy, min_curve.xy, max_curve.xy);
	fprintf( pldata,"strain zx            %5d     %5d   %14.6e %14.6e\n", min_strain_node.zx,
		max_strain_node.zx, min_strain.zx, max_strain.zx);
	fprintf( pldata,"strain yz            %5d     %5d   %14.6e %14.6e\n", min_strain_node.yz,
		max_strain_node.yz, min_strain.yz, max_strain.yz);
	fprintf( pldata,"curve I              %5d     %5d   %14.6e %14.6e\n", min_curve_node.I,
		max_curve_node.I, min_curve.I, max_curve.I);
	fprintf( pldata,"curve II             %5d     %5d   %14.6e %14.6e\n", min_curve_node.II,
		max_curve_node.II, min_curve.II, max_curve.II);
	fprintf( pldata,"\n");
	fprintf( pldata,"Orthographic viewport parameters(right, left, top, bootom, near, far)\n ");
	fprintf( pldata,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ortho_right, ortho_left,
		ortho_top, ortho_bottom, near, 1000.0);
	fprintf( pldata,"Perspective viewport parameters( mesh width and height)\n ");
	fprintf( pldata,"%6d %6d\n", mesh_width, mesh_height);
	fprintf( pldata,"Step sizes in x, y, z\n ");
	fprintf( pldata,"%14.6e %14.6e %14.6e\n",step_sizex, step_sizey, step_sizez);
	fprintf( pldata,"Amplification size\n ");
	fprintf( pldata,"%14.6e\n",amplify_step0);

	fclose( pldata );

  	return 1;    /* ANSI C requires main to return int. */
}

