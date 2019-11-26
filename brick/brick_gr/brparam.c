/*
    This program calculates and writes the parameters for
    the FEM GUI for brick elements.
  
   			Last Update 1/23/02

    SLFFEA source file
    Version:  1.3
    Copyright (C) 1999, 2000, 2001, 2002  San Le 

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
#include "../brick/brconst.h"
#include "../brick/brstruct.h"
#include "brstrcgr.h"
#include "../../common_gr/control.h"

#define init_far0      -2.0

extern int nmat, numnp, numel, dof;
extern double step_sizex, step_sizey, step_sizez;
extern double left, right, top, bottom, near, far, fscale, coord_rescale;
extern double cross_sec_left_right, cross_sec_up_down, cross_sec_in_out,
        cross_sec_left_right0, cross_sec_up_down0, cross_sec_in_out0;
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
	max_Uz, min_Uz, del_Uz, absolute_max_U, absolute_max_coord;


int brparameter( double *coord, SDIM *strain_node, SDIM *stress_node, double *U)
{
        int i, j, check;
	int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min,
		node_Uz_max, node_Uz_min; 
 	ISDIM max_stress_node, min_stress_node, max_strain_node, min_strain_node;
	FILE *brdata;

/*   brdata contains all the parameters and extreme values 
*/
	brdata = fopen( "brview.dat","w" );

/* Initialize parameters */

	step_sizex = .1; step_sizey = .1; step_sizez = .1;
        init_right = - BIG; init_left = BIG;
	init_top = - BIG; init_bottom = BIG;
	init_near = - BIG; init_far = init_far0; true_far = BIG;
	max_Ux = - BIG; min_Ux = BIG;
	max_Uy = - BIG; min_Uy = BIG;
	max_Uz = - BIG; min_Uz = BIG;

	node_Ux_max = 0; node_Ux_min = IBIG;
	node_Uy_max = 0; node_Uy_min = IBIG;
	node_Uz_max = 0; node_Uz_min = IBIG;
	 
 	max_strain_node.xx = 0; min_strain_node.xx = IBIG;
 	max_strain_node.yy = 0; min_strain_node.yy = IBIG;
 	max_strain_node.zz = 0; min_strain_node.zz = IBIG;
 	max_strain_node.xy = 0; min_strain_node.xy = IBIG;
 	max_strain_node.zx = 0; min_strain_node.zx = IBIG;
 	max_strain_node.yz = 0; min_strain_node.yz = IBIG;
 	max_strain_node.I = 0; min_strain_node.I = IBIG;
 	max_strain_node.II = 0; min_strain_node.II = IBIG;
 	max_strain_node.III = 0; min_strain_node.III = IBIG;

 	max_stress_node.xx = 0; min_stress_node.xx = IBIG;
 	max_stress_node.yy = 0; min_stress_node.yy = IBIG;
 	max_stress_node.zz = 0; min_stress_node.zz = IBIG;
 	max_stress_node.xy = 0; min_stress_node.xy = IBIG;
 	max_stress_node.zx = 0; min_stress_node.zx = IBIG;
 	max_stress_node.yz = 0; min_stress_node.yz = IBIG;
 	max_stress_node.I = 0; min_stress_node.I = IBIG;
 	max_stress_node.II = 0; min_stress_node.II = IBIG;
 	max_stress_node.III = 0; min_stress_node.III = IBIG;

/* Initialize largest and smallest strains */

        max_strain.xx = - BIG; min_strain.xx = BIG;
        max_strain.yy = - BIG; min_strain.yy = BIG;
        max_strain.zz = - BIG; min_strain.zz = BIG;
        max_strain.xy = - BIG; min_strain.xy = BIG;
        max_strain.zx = - BIG; min_strain.zx = BIG;
        max_strain.yz = - BIG; min_strain.yz = BIG;
        max_strain.I = - BIG; min_strain.I = BIG;
        max_strain.II = - BIG; min_strain.II = BIG;
        max_strain.III = - BIG; min_strain.III = BIG;

/* Initialize largest and smallest stresses */

        max_stress.xx = - BIG; min_stress.xx = BIG;
        max_stress.yy = - BIG; min_stress.yy = BIG;
        max_stress.zz = - BIG; min_stress.zz = BIG;
        max_stress.xy = - BIG; min_stress.xy = BIG;
        max_stress.zx = - BIG; min_stress.zx = BIG;
        max_stress.yz = - BIG; min_stress.yz = BIG;
        max_stress.I = - BIG; min_stress.I = BIG;
        max_stress.II = - BIG; min_stress.II = BIG;
        max_stress.III = - BIG; min_stress.III = BIG;

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

                if( init_near < *(coord+nsd*i+2))
                        init_near=*(coord+nsd*i+2);

                if( true_far > *(coord+nsd*i+2))
			true_far=*(coord+nsd*i+2);

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
                if( max_Uz < *(U+ndof*i+2))
		{
                        max_Uz=*(U+ndof*i+2);
			node_Uz_max = i;
		}
                if( min_Uz > *(U+ndof*i+2))
		{
                        min_Uz=*(U+ndof*i+2);
			node_Uz_min = i;
		}
        }


/* Because Mesa has problems with Meshes that have dimensions larger than 1000
   or smaller than .1, I am rescaling everything so that things are on the order
   of 10.0.
*/

	absolute_max_coord = fabs(init_left);
	if(absolute_max_coord < fabs(init_right)) absolute_max_coord = fabs(init_right);
	if(absolute_max_coord < fabs(init_bottom)) absolute_max_coord = fabs(init_bottom);
	if(absolute_max_coord < fabs(init_top)) absolute_max_coord = fabs(init_top);
	if(absolute_max_coord < fabs(init_near)) absolute_max_coord = fabs(init_near);
	if(absolute_max_coord < fabs(true_far)) absolute_max_coord = fabs(true_far);

	coord_rescale = 1.0;
	if( absolute_max_coord > 10.0 )     coord_rescale = 10.0;
	if( absolute_max_coord > 100.0 )    coord_rescale = 100.0;
	if( absolute_max_coord > 1000.0 )   coord_rescale = 1000.0;
	if( absolute_max_coord > 10000.0 )  coord_rescale = 10000.0;
	if( absolute_max_coord > 100000.0 ) coord_rescale = 100000.0;

	if( absolute_max_coord < 1.0 )     coord_rescale = 0.1;
	if( absolute_max_coord < 0.1 )     coord_rescale = 0.01;
	if( absolute_max_coord < 0.01 )    coord_rescale = 0.001;
	if( absolute_max_coord < 0.001 )   coord_rescale = 0.0001;
	if( absolute_max_coord < 0.0001 )  coord_rescale = 0.00001;
	if( absolute_max_coord < 0.00001 ) coord_rescale = 0.000001;

/* Rescale coordinates and displacements */

	if( coord_rescale > 1.01 || coord_rescale < .99 )
	{
		for( i = 0; i < numnp; ++i )
		{
                        *(coord+nsd*i) /= coord_rescale;
                        *(coord+nsd*i+1) /= coord_rescale;
                        *(coord+nsd*i+2) /= coord_rescale;

                        *(U+ndof*i) /= coord_rescale;
                        *(U+ndof*i+1) /= coord_rescale;
                        *(U+ndof*i+2) /= coord_rescale;

		}

		init_left /= coord_rescale;
		init_right /= coord_rescale;
		init_bottom /= coord_rescale;
		init_top /= coord_rescale;
		init_near /= coord_rescale;
		true_far /= coord_rescale;

		min_Ux /= coord_rescale;
		max_Ux /= coord_rescale;
		min_Uy /= coord_rescale;
		max_Uy /= coord_rescale;
		min_Uz /= coord_rescale;
		max_Uz /= coord_rescale;
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
                if( max_strain.zz < strain_node[i].zz )
		{
                        max_strain.zz = strain_node[i].zz;
			max_strain_node.zz = i;
		}
                if( min_strain.zz > strain_node[i].zz )
		{
                        min_strain.zz = strain_node[i].zz;
			min_strain_node.zz = i;
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
                if( max_strain.III < strain_node[i].III )
		{
                        max_strain.III = strain_node[i].III;
			max_strain_node.III = i;
		}
                if( min_strain.III > strain_node[i].III )
		{
                        min_strain.III = strain_node[i].III;
			min_strain_node.III = i;
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
                if( max_stress.zz < stress_node[i].zz )
		{
                        max_stress.zz = stress_node[i].zz;
			max_stress_node.zz = i;
		}
                if( min_stress.zz > stress_node[i].zz )
		{
                        min_stress.zz = stress_node[i].zz;
			min_stress_node.zz = i;
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
                if( max_stress.III < stress_node[i].III )
		{
                        max_stress.III = stress_node[i].III;
			max_stress_node.III = i;
		}
                if( min_stress.III > stress_node[i].III )
		{
                        min_stress.III = stress_node[i].III;
			min_stress_node.III = i;
		}
        }



/* Set the axes parameters */

	AxisMax_x = 1.2*init_right + 1.0;
	AxisMax_y = 1.2*init_top + 1.0;
	AxisMax_z = 1.2*init_near + 1.0;
	AxisMin_x = 1.2*init_left - 1.0;
	AxisMin_y = 1.2*init_bottom - 1.0;
	AxisMin_z = 1.2*true_far - 1.0;

	if( AxisMax_x < 0.0 )
		AxisMax_x = -.1*AxisMin_x;
	if( AxisMax_y < 0.0 )
		AxisMax_y = -.1*AxisMin_y;
	if( AxisMax_z < 0.0 )
		AxisMax_z = -.1*AxisMin_z;
	if( AxisMin_x > 0.0 )
		AxisMin_x = -.1*AxisMax_x;
	if( AxisMin_y > 0.0 )
		AxisMin_y = -.1*AxisMax_y;
	if( AxisMin_z > 0.0 )
		AxisMin_z = -.1*AxisMax_z;

	IAxisMin_x = (double)(int)AxisMin_x; 
	IAxisMin_y = (double)(int)AxisMin_y;
	IAxisMin_z = (double)(int)AxisMin_z;

	AxisLength_x = (int)(AxisMax_x - AxisMin_x);
	AxisLength_y = (int)(AxisMax_y - AxisMin_y);
	AxisLength_z = (int)(AxisMax_z - AxisMin_z);
	AxisLength_max = AxisLength_x;
	if( AxisLength_max < AxisLength_y )
		AxisLength_max = AxisLength_y;
	if( AxisLength_max < AxisLength_z )
		AxisLength_max = AxisLength_z;

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
	step_sizez = (near - far) / 5.0;

	left_right0 = -(left + right) / 2.0;
	up_down0 = - (top + bottom ) / 2.0;
	/*in_out  = (far + near ) / 2.0 - 5.0;*/
	in_out0 = far - 20.0*AxisLength_max;
	left_right = left_right0;
	up_down = up_down0;
	in_out = in_out0;

/* Set the Cross Section Plane parameters */

	cross_sec_left_right0 = AxisMax_x;
	cross_sec_up_down0 = AxisMax_y;
	cross_sec_in_out0 = AxisMax_z;

	cross_sec_left_right = cross_sec_left_right0;
	cross_sec_up_down = cross_sec_up_down0;
	cross_sec_in_out = cross_sec_in_out0;

    	mesh_width = mesh_width0;
    	mesh_height = mesh_height0;

/* Print the above data in the file "brview.dat" */

	fprintf( brdata, "                            node\n");
	fprintf( brdata, "                          min  max       min            max\n");
	fprintf( brdata,"displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux*coord_rescale, max_Ux*coord_rescale);
	fprintf( brdata,"displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy*coord_rescale, max_Uy*coord_rescale);
	fprintf( brdata,"displacement Uz        %5d %5d   %14.6e %14.6e\n", node_Uz_min,
		node_Uz_max, min_Uz*coord_rescale, max_Uz*coord_rescale);
	fprintf( brdata,"\n");
	fprintf( brdata, "                            node\n");
	fprintf( brdata, "                        min       max         min           max\n");
	fprintf( brdata,"stress xx            %5d     %5d   %14.6e %14.6e\n", min_stress_node.xx,
		max_stress_node.xx, min_stress.xx, max_stress.xx);
	fprintf( brdata,"stress yy            %5d     %5d   %14.6e %14.6e\n", min_stress_node.yy,
		max_stress_node.yy, min_stress.yy, max_stress.yy);
	fprintf( brdata,"stress zz            %5d     %5d   %14.6e %14.6e\n", min_stress_node.zz,
		max_stress_node.zz, min_stress.zz, max_stress.zz);
	fprintf( brdata,"stress xy            %5d     %5d   %14.6e %14.6e\n", min_stress_node.xy,
		max_stress_node.xy, min_stress.xy, max_stress.xy);
	fprintf( brdata,"stress zx            %5d     %5d   %14.6e %14.6e\n", min_stress_node.zx,
		max_stress_node.zx, min_stress.zx, max_stress.zx);
	fprintf( brdata,"stress yz            %5d     %5d   %14.6e %14.6e\n", min_stress_node.yz,
		max_stress_node.yz, min_stress.yz, max_stress.yz);
	fprintf( brdata,"stress I             %5d     %5d   %14.6e %14.6e\n", min_stress_node.I,
		max_stress_node.I, min_stress.I, max_stress.I);
	fprintf( brdata,"stress II            %5d     %5d   %14.6e %14.6e\n", min_stress_node.II,
		max_stress_node.II, min_stress.II, max_stress.II);
	fprintf( brdata,"stress III           %5d     %5d   %14.6e %14.6e\n", min_stress_node.III,
		max_stress_node.III, min_stress.III, max_stress.III);
	fprintf( brdata,"\n");
	fprintf( brdata,"strain xx            %5d     %5d   %14.6e %14.6e\n", min_strain_node.xx,
		max_strain_node.xx, min_strain.xx, max_strain.xx);
	fprintf( brdata,"strain yy            %5d     %5d   %14.6e %14.6e\n", min_strain_node.yy,
		max_strain_node.yy, min_strain.yy, max_strain.yy);
	fprintf( brdata,"strain zz            %5d     %5d   %14.6e %14.6e\n", min_strain_node.zz,
		max_strain_node.zz, min_strain.zz, max_strain.zz);
	fprintf( brdata,"strain xy            %5d     %5d   %14.6e %14.6e\n", min_strain_node.xy,
		max_strain_node.xy, min_strain.xy, max_strain.xy);
	fprintf( brdata,"strain zx            %5d     %5d   %14.6e %14.6e\n", min_strain_node.zx,
		max_strain_node.zx, min_strain.zx, max_strain.zx);
	fprintf( brdata,"strain yz            %5d     %5d   %14.6e %14.6e\n", min_strain_node.yz,
		max_strain_node.yz, min_strain.yz, max_strain.yz);
	fprintf( brdata,"strain I             %5d     %5d   %14.6e %14.6e\n", min_strain_node.I,
		max_strain_node.I, min_strain.I, max_strain.I);
	fprintf( brdata,"strain II            %5d     %5d   %14.6e %14.6e\n", min_strain_node.II,
		max_strain_node.II, min_strain.II, max_strain.II);
	fprintf( brdata,"strain III           %5d     %5d   %14.6e %14.6e\n", min_strain_node.III,
		max_strain_node.III, min_strain.III, max_strain.III);
	fprintf( brdata,"\n");
	fprintf( brdata,"Orthographic viewport parameters(right, left, top, bootom, near, far)\n ");
	fprintf( brdata,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ortho_right, ortho_left,
		ortho_top, ortho_bottom, near, 1000.0);
	fprintf( brdata,"Perspective viewport parameters( mesh width and height)\n ");
	fprintf( brdata,"%6d %6d\n", mesh_width, mesh_height);
	fprintf( brdata,"Step sizes in x, y, z\n ");
	fprintf( brdata,"%14.6e %14.6e %14.6e\n",step_sizex, step_sizey, step_sizez);
	fprintf( brdata,"Amplification size\n ");
	fprintf( brdata,"%14.6e\n",amplify_step0);

	fclose( brdata );

  	return 1;    /* ANSI C requires main to return int. */
}
