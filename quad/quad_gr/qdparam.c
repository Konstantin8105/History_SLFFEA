/*
    This program calculates and writes the parameters for
    the FEM GUI for quad elements.
  
   			Last Update 6/4/00

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
#include "../quad/qdconst.h"
#include "../quad/qdstruct.h"
#include "qdstrcgr.h"
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
extern SDIM del_stress, del_strain, max_stress, min_stress,
	max_strain, min_strain;
extern double max_Ux, min_Ux, del_Ux, max_Uy, min_Uy, del_Uy,
	max_Uz, min_Uz, del_Uz, absolute_max_U;

int qdparameter( double *coord, SDIM *strain_node, SDIM *stress_node, double *U )
{
        int i, j, check;
	int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min;
	ISDIM max_stress_node, min_stress_node, max_strain_node, min_strain_node;
	FILE *qddata;
	
/*   qddata contains all the parameters and extreme values
*/
	qddata = fopen( "qdview.dat","w" );

/* Initialize parameters */

	step_sizex = .1; step_sizey = .1; step_sizez = .1;
	init_right = - BIG; init_left = BIG;
	init_top = - BIG; init_bottom = BIG;
	init_near = 0.0; init_far = 0.0; true_far = 0.0;
	max_Ux = - BIG; min_Ux = BIG;
	max_Uy = - BIG; min_Uy = BIG;
	max_Uz = 0.0; min_Uz = 0.0;

	node_Ux_max = 0;
	node_Ux_min = IBIG;
	node_Uy_max = 0;
	node_Uy_min = IBIG;

	max_strain_node.xx = 0; min_strain_node.xx = IBIG;
	max_strain_node.yy = 0; min_strain_node.yy = IBIG;
	max_strain_node.xy = 0; min_strain_node.xy = IBIG;
	max_strain_node.I = 0; min_strain_node.I = IBIG;
	max_strain_node.II = 0; min_strain_node.II = IBIG;

	max_stress_node.xx = 0; min_stress_node.xx = IBIG;
	max_stress_node.yy = 0; min_stress_node.yy = IBIG;
	max_stress_node.xy = 0; min_stress_node.xy = IBIG;
	max_stress_node.I = 0; min_stress_node.I = IBIG;
	max_stress_node.II = 0; min_stress_node.II = IBIG;

/* Initialize largest and smallest strains */

	max_strain.xx = - BIG; min_strain.xx = BIG;
	max_strain.yy = - BIG; min_strain.yy = BIG;
	max_strain.xy = - BIG; min_strain.xy = BIG;
	max_strain.I = - BIG; min_strain.I = BIG;
	max_strain.II = - BIG; min_strain.II = BIG;

/* Initialize largest and smallest stresses */

	max_stress.xx = - BIG; min_stress.xx = BIG;
	max_stress.yy = - BIG; min_stress.yy = BIG;
	max_stress.xy = - BIG; min_stress.xy = BIG;
	max_stress.I = - BIG; min_stress.I = BIG;
	max_stress.II = - BIG; min_stress.II = BIG;

/* Search for extreme values */
 
/* Search for extreme values of displacement and nodal point */

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

/* Search for extreme nodal displacements */

		if( max_Ux < *(U+ndof*i))
		{
			max_Ux=*(U+ndof*i);
			node_Ux_max = i;
		}
		if( min_Ux > *(U+ndof*i))
		{
			min_Ux=*(U+ndof*i);
			node_Ux_min = i;
		}
		if( max_Uy < *(U+ndof*i+1))
		{
			max_Uy=*(U+ndof*i+1);
			node_Uy_max = i;
		}
		if( min_Uy > *(U+ndof*i+1))
		{
			min_Uy=*(U+ndof*i+1);
			node_Uy_min = i;
		}
        }

/* Search for largest absolute value U */

	absolute_max_U = fabs(min_Ux);
        if( absolute_max_U < fabs(max_Ux))
		absolute_max_U = fabs(max_Ux);
        if( absolute_max_U < fabs(min_Uy))
		absolute_max_U = fabs(min_Uy);
        if( absolute_max_U < fabs(max_Uy))
		absolute_max_U = fabs(max_Uy);
        if( absolute_max_U < fabs(min_Uz))
		absolute_max_U = fabs(min_Uz);
        if( absolute_max_U < fabs(max_Uz))
		absolute_max_U = fabs(max_Uz);

        if( init_far > true_far)
		init_far=true_far;

        for( i = 0; i < numnp; ++i )
        {

/* Find extreme strains */

                if( max_strain.xx < strain_node[i].xx )
		{
                        max_strain.xx = strain_node[i].xx;
			max_strain_node.xx = i;
		}
                if( min_strain.xx > strain_node[i].xx )
		{
                        min_strain.xx = strain_node[i].xx;
			min_strain_node.xx = i;
		}
                if( max_strain.yy < strain_node[i].yy )
		{
                        max_strain.yy = strain_node[i].yy;
			max_strain_node.yy = i;
		}
                if( min_strain.yy > strain_node[i].yy )
		{
                        min_strain.yy = strain_node[i].yy;
			min_strain_node.yy = i;
		}
                if( max_strain.xy < strain_node[i].xy )
		{
                        max_strain.xy = strain_node[i].xy;
			max_strain_node.xy = i;
		}
                if( min_strain.xy > strain_node[i].xy )
		{
                        min_strain.xy = strain_node[i].xy;
			min_strain_node.xy = i;
		}
                if( max_strain.I < strain_node[i].I )
		{
                        max_strain.I = strain_node[i].I;
			max_strain_node.I = i;
		}
                if( min_strain.I > strain_node[i].I )
		{
                        min_strain.I = strain_node[i].I;
			min_strain_node.I = i;
		}
                if( max_strain.II < strain_node[i].II )
		{
                        max_strain.II = strain_node[i].II;
			max_strain_node.II = i;
		}
                if( min_strain.II > strain_node[i].II )
		{
                        min_strain.II = strain_node[i].II;
			min_strain_node.II = i;
		}
/* Find extreme stresses */

                if( max_stress.xx < stress_node[i].xx )
		{
                        max_stress.xx = stress_node[i].xx;
			max_stress_node.xx = i;
		}
                if( min_stress.xx > stress_node[i].xx )
		{
                        min_stress.xx = stress_node[i].xx;
			min_stress_node.xx = i;
		}
                if( max_stress.yy < stress_node[i].yy )
		{
                        max_stress.yy = stress_node[i].yy;
			max_stress_node.yy = i;
		}
                if( min_stress.yy > stress_node[i].yy )
		{
                        min_stress.yy = stress_node[i].yy;
			min_stress_node.yy = i;
		}
                if( max_stress.xy < stress_node[i].xy )
		{
                        max_stress.xy = stress_node[i].xy;
			max_stress_node.xy = i;
		}
                if( min_stress.xy > stress_node[i].xy )
		{
                        min_stress.xy = stress_node[i].xy;
			min_stress_node.xy = i;
		}
                if( max_stress.I < stress_node[i].I )
		{
                        max_stress.I = stress_node[i].I;
			max_stress_node.I = i;
		}
                if( min_stress.I > stress_node[i].I )
		{
                        min_stress.I = stress_node[i].I;
			min_stress_node.I = i;
		}
                if( max_stress.II < stress_node[i].II )
		{
                        max_stress.II = stress_node[i].II;
			max_stress_node.II = i;
		}
                if( min_stress.II > stress_node[i].II )
		{
                        min_stress.II = stress_node[i].II;
			min_stress_node.II = i;
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
	/*printf(" AxisPoint_step %10.5e\n",AxisPoint_step);*/

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

/* Print the above data in the file "qdview.dat" */

	fprintf( qddata, "                            node\n");
	fprintf( qddata, "                          min  max       min            max\n");
	fprintf( qddata,"displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux, max_Ux);
	fprintf( qddata,"displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy, max_Uy);
	fprintf( qddata,"\n");
	fprintf( qddata, "                            node\n");
	fprintf( qddata, "                        min       max         min           max\n");
	fprintf( qddata,"stress xx            %5d     %5d   %14.6e %14.6e\n", min_stress_node.xx,
		max_stress_node.xx, min_stress.xx, max_stress.xx);
	fprintf( qddata,"stress yy            %5d     %5d   %14.6e %14.6e\n", min_stress_node.yy,
		max_stress_node.yy, min_stress.yy, max_stress.yy);
	fprintf( qddata,"stress xy            %5d     %5d   %14.6e %14.6e\n", min_stress_node.xy,
		max_stress_node.xy, min_stress.xy, max_stress.xy);
	fprintf( qddata,"stress I             %5d     %5d   %14.6e %14.6e\n", min_stress_node.I,
		max_stress_node.I, min_stress.I, max_stress.I);
	fprintf( qddata,"stress II            %5d     %5d   %14.6e %14.6e\n", min_stress_node.II,
		max_stress_node.II, min_stress.II, max_stress.II);
	fprintf( qddata,"\n");
	fprintf( qddata,"strain xx            %5d     %5d   %14.6e %14.6e\n", min_strain_node.xx,
		max_strain_node.xx, min_strain.xx, max_strain.xx);
	fprintf( qddata,"strain yy            %5d     %5d   %14.6e %14.6e\n", min_strain_node.yy,
		max_strain_node.yy, min_strain.yy, max_strain.yy);
	fprintf( qddata,"strain xy            %5d     %5d   %14.6e %14.6e\n", min_strain_node.xy,
		max_strain_node.xy, min_strain.xy, max_strain.xy);
	fprintf( qddata,"strain I             %5d     %5d   %14.6e %14.6e\n", min_strain_node.I,
		max_strain_node.I, min_strain.I, max_strain.I);
	fprintf( qddata,"strain II            %5d     %5d   %14.6e %14.6e\n", min_strain_node.II,
		max_strain_node.II, min_strain.II, max_strain.II);
	fprintf( qddata,"\n");
	fprintf( qddata,"Orthographic viewport parameters(right, left, top, bootom, near, far)\n ");
	fprintf( qddata,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ortho_right, ortho_left,
		ortho_top, ortho_bottom, near, 1000.0);
	fprintf( qddata,"Perspective viewport parameters( mesh width and height)\n ");
	fprintf( qddata,"%6d %6d\n", mesh_width, mesh_height);
	fprintf( qddata,"Step sizes in x, y, z\n ");
	fprintf( qddata,"%14.6e %14.6e %14.6e\n",step_sizex, step_sizey, step_sizez);
	fprintf( qddata,"Amplification size\n ");
	fprintf( qddata,"%14.6e\n",amplify_step0);

	fclose( qddata );

  	return 1;    /* ANSI C requires main to return int. */
}
