/*
    This program calculates and writes the parameters for
    the FEM GUI for truss elements.
  
   			Last Update 6/12/01

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001  San Le 

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
#include "../truss/tsconst.h"
#include "../truss/tsstruct.h"
#include "tsstrcgr.h"
#include "../../common_gr/control.h"

#define init_far0      -2.0

extern int nmat, numnp, numel, dof;
extern double step_sizex, step_sizey, step_sizez;
extern double left, right, top, bottom, near, far, fscale;
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
extern STRESS del_stress, max_stress, min_stress;
extern STRAIN del_strain, max_strain, min_strain;
extern double max_Ux, min_Ux, del_Ux, max_Uy, min_Uy, del_Uy,
	max_Uz, min_Uz, del_Uz, absolute_max_U;

int tsparameter( double *coord, STRAIN *strain, STRESS *stress, double *U )
{
        int i, j, check;
        int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min, node_Uz_max, node_Uz_min;
        ISDIM max_stress_el, min_stress_el, max_stress_integ, min_stress_integ,
                max_strain_el, min_strain_el, max_strain_integ, min_strain_integ;
        FILE *tsdata;

/*   shdata contains all the parameters and extreme values
*/
        tsdata = fopen( "tsview.dat","w" );

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

	max_strain_el.xx = 0; min_strain_el.xx = IBIG;
	max_strain_integ.xx = 0; min_strain_integ.xx = IBIG;
	max_stress_el.xx = 0; min_stress_el.xx = IBIG;
	max_stress_integ.xx = 0; min_stress_integ.xx = IBIG;

/* Initialize for largest and smallest strains */

        max_strain.xx = - BIG; min_strain.xx = BIG;

/* Initialize largest and smallest stresses */

        max_stress.xx = - BIG; min_stress.xx = BIG;

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

/* Search for largest absolute value of displacement U */

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

	max_strain_integ.xx = 0;
	min_strain_integ.xx = 0;
	max_stress_integ.xx = 0;
	min_stress_integ.xx = 0;

	for( i = 0; i < numel; ++i )
	{

/* Find extreme strains */

		if( max_strain.xx < strain[i].xx )
		{
			max_strain.xx = strain[i].xx;
			max_strain_el.xx = i;
		}
		if( min_strain.xx > strain[i].xx )
		{
			min_strain.xx = strain[i].xx;
			min_strain_el.xx = i;
		}
/* Find extreme stresses */

		if( max_stress.xx < stress[i].xx )
		{
			max_stress.xx = stress[i].xx;
			max_stress_el.xx = i;
		}
		if( min_stress.xx > stress[i].xx )
		{
			min_stress.xx = stress[i].xx;
			min_stress_el.xx = i;
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
	/*printf(" amplify_step0 %10.5e %10.5e %10.5e\n",amplify_step0,
		AxisLength_max,absolute_max_U);*/

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

    	mesh_width=mesh_width0;
    	mesh_height=mesh_height0;


/* Print the above data in the file "tsview.dat" */

	fprintf( tsdata, "                            node\n");
	fprintf( tsdata, "                          min  max       min            max\n");
	fprintf( tsdata,"displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux, max_Ux);
	fprintf( tsdata,"displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy, max_Uy);
	fprintf( tsdata,"displacement Uz        %5d %5d   %14.6e %14.6e\n", node_Uz_min,
		node_Uz_max, min_Uz, max_Uz);
	fprintf( tsdata,"\n");
	fprintf( tsdata, "                        el. gauss pt.\n");
	fprintf( tsdata, "                        min       max         min           max\n");
	fprintf( tsdata,"stress xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_stress_el.xx,
		min_stress_integ.xx, max_stress_el.xx, max_stress_integ.xx,
		min_stress.xx, max_stress.xx);
	fprintf( tsdata,"\n");
	fprintf( tsdata,"strain xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_strain_el.xx,
		min_strain_integ.xx, max_strain_el.xx, max_strain_integ.xx,
		min_strain.xx, max_strain.xx);
	fprintf( tsdata,"\n");
	fprintf( tsdata,"Orthographic viewport parameters(right, left, top, bootom, near, far)\n ");
	fprintf( tsdata,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ortho_right, ortho_left,
		ortho_top, ortho_bottom, near, 1000.0);
	fprintf( tsdata,"Perspective viewport parameters( mesh width and height)\n ");
	fprintf( tsdata,"%6d %6d\n", mesh_width, mesh_height);
	fprintf( tsdata,"Step sizes in x, y, z\n ");
	fprintf( tsdata,"%14.6e %14.6e %14.6e\n",step_sizex, step_sizey, step_sizez);
	fprintf( tsdata,"Amplification size\n ");
	fprintf( tsdata,"%14.6e\n",amplify_step0);

	fclose( tsdata );

	return 1;    /* ANSI C requires main to return int. */
}

