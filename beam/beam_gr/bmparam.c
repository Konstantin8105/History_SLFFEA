/*
    This program calculates and writes the parameters for
    the FEM GUI for beam elements.
  
   			Last Update 3/17/00

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
#include "../beam/bmconst.h"
#include "../beam/bmstruct.h"
#include "bmstrcgr.h"
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
	max_Uphi_z, min_Uphi_z, del_Uphi_z,
	max_Ux, min_Ux, del_Ux, max_Uy, min_Uy, del_Uy,
	max_Uz, min_Uz, del_Uz, absolute_max_U;

int bmparameter(double *coord, CURVATURE *curve, MOMENT *moment,
	STRAIN *strain, STRESS *stress, double *U )
{
        int i, j, check; 
	int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min, node_Uz_max, node_Uz_min,
		node_Uphi_x_max, node_Uphi_x_min, node_Uphi_y_max, node_Uphi_y_min,
		node_Uphi_z_max, node_Uphi_z_min;
	IMDIM max_moment_el, min_moment_el, max_moment_integ, min_moment_integ,
		max_curve_el, min_curve_el, max_curve_integ, min_curve_integ;
	ISDIM max_stress_el, min_stress_el, max_stress_integ, min_stress_integ,
		max_strain_el, min_strain_el, max_strain_integ, min_strain_integ;
	FILE *bmdata;

/*   bmdata contains all the parameters and extreme values
*/
        bmdata = fopen( "bmview.dat","w" );

/* Initialize parameters */

	step_sizex = .1; step_sizey = .1; step_sizez = .1;
        init_right = - BIG; init_left = BIG;
	init_top = - BIG; init_bottom = BIG;
	init_near = - BIG; init_far = init_far0; true_far = BIG;
	max_Ux = - BIG; min_Ux = BIG;
	max_Uy = - BIG; min_Uy = BIG;
	max_Uz = - BIG; min_Uz = BIG;
	max_Uphi_x = - BIG; min_Uphi_x = BIG;
	max_Uphi_y = - BIG; min_Uphi_y = BIG;
	max_Uphi_z = - BIG; min_Uphi_z = BIG;

        node_Ux_max = 0; node_Ux_min = IBIG;
        node_Uy_max = 0; node_Uy_min = IBIG;
        node_Uz_max = 0; node_Uz_min = IBIG;

        node_Uphi_x_max = 0; node_Uphi_x_min = IBIG;
        node_Uphi_y_max = 0; node_Uphi_y_min = IBIG;
        node_Uphi_z_max = 0; node_Uphi_z_min = IBIG;

        max_curve_el.xx = 0; min_curve_el.xx = IBIG;
        max_curve_el.yy = 0; min_curve_el.yy = IBIG;
        max_curve_el.zz = 0; min_curve_el.zz = IBIG;
        max_strain_el.xx = 0; min_strain_el.xx = IBIG;
        max_curve_integ.xx = 0; min_curve_integ.xx = IBIG;
        max_curve_integ.yy = 0; min_curve_integ.yy = IBIG;
        max_curve_integ.zz = 0; min_curve_integ.zz = IBIG;
        max_strain_integ.xx = 0; min_strain_integ.xx = IBIG;

        max_moment_el.xx = 0; min_moment_el.xx = IBIG;
        max_moment_el.yy = 0; min_moment_el.yy = IBIG;
        max_moment_el.zz = 0; min_moment_el.zz = IBIG;
        max_stress_el.xx = 0; min_stress_el.xx = IBIG;
        max_moment_integ.xx = 0; min_moment_integ.xx = IBIG;
        max_moment_integ.yy = 0; min_moment_integ.yy = IBIG;
        max_moment_integ.zz = 0; min_moment_integ.zz = IBIG;
        max_stress_integ.xx = 0; min_stress_integ.xx = IBIG;

/* Initialize for largest and smallest curvatures and strains */

        max_curve.xx = - BIG; min_curve.xx = BIG;
        max_curve.yy = - BIG; min_curve.yy = BIG;
        max_curve.zz = - BIG; min_curve.zz = BIG;
        max_strain.xx = - BIG; min_strain.xx = BIG;

/* Initialize largest and smallest moments and stresses */

        max_moment.xx = - BIG; min_moment.xx = BIG;
        max_moment.yy = - BIG; min_moment.yy = BIG;
        max_moment.zz = - BIG; min_moment.zz = BIG;
        max_stress.xx = - BIG; min_stress.xx = BIG;

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

		if( init_near < *(coord+nsd*i+2))
			init_near=*(coord+nsd*i+2);

		if( true_far > *(coord+nsd*i+2))
			true_far=*(coord+nsd*i+2);

/* Search for extreme nodal displacements and angles */

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
                if( max_Uphi_x < *(U+ndof*i+3))
		{
                        max_Uphi_x=*(U+ndof*i+3);
			node_Uphi_x_max = i;
		}
                if( min_Uphi_x > *(U+ndof*i+3))
		{
                        min_Uphi_x=*(U+ndof*i+3);
			node_Uphi_x_min = i;
		}
                if( max_Uphi_y < *(U+ndof*i+4))
		{
                        max_Uphi_y=*(U+ndof*i+4);
			node_Uphi_y_max = i;
		}
                if( min_Uphi_y > *(U+ndof*i+4))
		{
                        min_Uphi_y=*(U+ndof*i+4);
			node_Uphi_y_min = i;
		}
                if( max_Uphi_z < *(U+ndof*i+5))
		{
                        max_Uphi_z=*(U+ndof*i+5);
			node_Uphi_z_max = i;
		}
                if( min_Uphi_z > *(U+ndof*i+5))
		{
                        min_Uphi_z=*(U+ndof*i+5);
			node_Uphi_z_min = i;
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

        for( i = 0; i < numel; ++i )
        {
            for( j = 0; j < num_int; ++j )
            {

/* Find extreme curvatures and strains */

                if( max_curve.xx < curve[i].pt[j].xx )
		{
                        max_curve.xx = curve[i].pt[j].xx;
			max_curve_el.xx = i;
			max_curve_integ.xx = j;
		}
                if( min_curve.xx > curve[i].pt[j].xx )
		{
                        min_curve.xx = curve[i].pt[j].xx;
			min_curve_el.xx = i;
			min_curve_integ.xx = j;
		}
                if( max_curve.yy < curve[i].pt[j].yy )
		{
                        max_curve.yy = curve[i].pt[j].yy;
			max_curve_el.yy = i;
			max_curve_integ.yy = j;
		}
                if( min_curve.yy > curve[i].pt[j].yy )
		{
                        min_curve.yy = curve[i].pt[j].yy;
			min_curve_el.yy = i;
			min_curve_integ.yy = j;
		}
                if( max_curve.zz < curve[i].pt[j].zz )
		{
                        max_curve.zz = curve[i].pt[j].zz;
			max_curve_el.zz = i;
			max_curve_integ.zz = j;
		}
                if( min_curve.zz > curve[i].pt[j].zz )
		{
                        min_curve.zz = curve[i].pt[j].zz;
			min_curve_el.zz = i;
			min_curve_integ.zz = j;
		}
                if( max_strain.xx < strain[i].pt[j].xx )
		{
                        max_strain.xx = strain[i].pt[j].xx;
			max_strain_el.xx = i;
			max_strain_integ.xx = j;
		}
                if( min_strain.xx > strain[i].pt[j].xx )
		{
                        min_strain.xx = strain[i].pt[j].xx;
			min_strain_el.xx = i;
			min_strain_integ.xx = j;
		}
/* Find extreme moments and stresses */

                if( max_moment.xx < moment[i].pt[j].xx )
		{
                        max_moment.xx = moment[i].pt[j].xx;
			max_moment_el.xx = i;
			max_moment_integ.xx = j;
		}
                if( min_moment.xx > moment[i].pt[j].xx )
		{
                        min_moment.xx = moment[i].pt[j].xx;
			min_moment_el.xx = i;
			min_moment_integ.xx = j;
		}
                if( max_moment.yy < moment[i].pt[j].yy )
		{
                        max_moment.yy = moment[i].pt[j].yy;
			max_moment_el.yy = i;
			max_moment_integ.yy = j;
		}
                if( min_moment.yy > moment[i].pt[j].yy )
		{
                        min_moment.yy = moment[i].pt[j].yy;
			min_moment_el.yy = i;
			min_moment_integ.yy = j;
		}
                if( max_moment.zz < moment[i].pt[j].zz )
		{
                        max_moment.zz = moment[i].pt[j].zz;
			max_moment_el.zz = i;
			max_moment_integ.zz = j;
		}
                if( min_moment.zz > moment[i].pt[j].zz )
		{
                        min_moment.zz = moment[i].pt[j].zz;
			min_moment_el.zz = i;
			min_moment_integ.zz = j;
		}
                if( max_stress.xx < stress[i].pt[j].xx )
		{
                        max_stress.xx = stress[i].pt[j].xx;
			max_stress_el.xx = i;
			max_stress_integ.xx = j;
		}
                if( min_stress.xx > stress[i].pt[j].xx )
		{
                        min_stress.xx = stress[i].pt[j].xx;
			min_stress_el.xx = i;
			min_stress_integ.xx = j;
		}
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

	step_sizex = (right - left) / 20.0 + AxisLength_max;
	step_sizey = (top - bottom) / 20.0 + AxisLength_max;
	step_sizez = (near - far) / 5.0 + AxisLength_max;

	left_right0 = -(left + right) / 2.0;
	up_down0 = - (top + bottom ) / 2.0;
	/*in_out  = (far + near ) / 2.0 - 5.0;*/
	in_out0 = far - 20.0*AxisLength_max;
	left_right = left_right0;
	up_down = up_down0;
	in_out = in_out0;

    	mesh_width=mesh_width0;
    	mesh_height=mesh_height0;

/* Print the above data in the file "bmview.dat" */

	fprintf( bmdata, "                            node\n");
	fprintf( bmdata, "                          min  max       min            max\n");
	fprintf( bmdata,"displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux, max_Ux);
	fprintf( bmdata,"displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy, max_Uy);
	fprintf( bmdata,"displacement Uz        %5d %5d   %14.6e %14.6e\n", node_Uz_min,
		node_Uz_max, min_Uz, max_Uz);
	fprintf( bmdata,"angle phi x            %5d %5d   %14.6e %14.6e\n", node_Uphi_x_min,
		node_Uphi_x_max, min_Uphi_x, max_Uphi_x);
	fprintf( bmdata,"angle phi y            %5d %5d   %14.6e %14.6e\n", node_Uphi_y_min,
		node_Uphi_y_max, min_Uphi_y, max_Uphi_y);
	fprintf( bmdata,"angle phi z            %5d %5d   %14.6e %14.6e\n", node_Uphi_z_min,
		node_Uphi_z_max, min_Uphi_z, max_Uphi_z);
	fprintf( bmdata,"\n");
	fprintf( bmdata, "                        el. gauss pt.\n");
	fprintf( bmdata, "                        min       max         min           max\n");
	fprintf( bmdata,"moment xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_moment_el.xx,
		min_moment_integ.xx, max_moment_el.xx, max_moment_integ.xx,
		min_moment.xx, max_moment.xx);
	fprintf( bmdata,"moment yy            %5d %2d %5d %2d  %14.6e %14.6e\n", min_moment_el.yy,
		min_moment_integ.yy, max_moment_el.yy, max_moment_integ.yy,
		min_moment.yy, max_moment.yy);
	fprintf( bmdata,"moment zz            %5d %2d %5d %2d  %14.6e %14.6e\n", min_moment_el.zz,
		min_moment_integ.zz, max_moment_el.zz, max_moment_integ.zz,
		min_moment.zz, max_moment.zz);
	fprintf( bmdata,"stress xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_stress_el.xx,
		min_stress_integ.xx, max_stress_el.xx, max_stress_integ.xx,
		min_stress.xx, max_stress.xx);
	fprintf( bmdata,"\n");
	fprintf( bmdata,"curve xx             %5d %2d %5d %2d  %14.6e %14.6e\n", min_curve_el.xx,
		min_curve_integ.xx, max_curve_el.xx, max_curve_integ.xx,
		min_curve.xx, max_curve.xx);
	fprintf( bmdata,"curve yy             %5d %2d %5d %2d  %14.6e %14.6e\n", min_curve_el.yy,
		min_curve_integ.yy, max_curve_el.yy, max_curve_integ.yy,
		min_curve.yy, max_curve.yy);
	fprintf( bmdata,"curve zz             %5d %2d %5d %2d  %14.6e %14.6e\n", min_curve_el.zz,
		min_curve_integ.zz, max_curve_el.zz, max_curve_integ.zz,
		min_curve.zz, max_curve.zz);
	fprintf( bmdata,"strain xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_strain_el.xx,
		min_strain_integ.xx, max_strain_el.xx, max_strain_integ.xx,
		min_strain.xx, max_strain.xx);
	fprintf( bmdata,"\n");
	fprintf( bmdata,"Orthographic viewport parameters(right, left, top, bootom, near, far)\n ");
	fprintf( bmdata,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ortho_right, ortho_left,
		ortho_top, ortho_bottom, near, 1000.0);
	fprintf( bmdata,"Perspective viewport parameters( mesh width and height)\n ");
	fprintf( bmdata,"%6d %6d\n", mesh_width, mesh_height);
	fprintf( bmdata,"Step sizes in x, y, z\n ");
	fprintf( bmdata,"%14.6e %14.6e %14.6e\n",step_sizex, step_sizey, step_sizez);
	fprintf( bmdata,"Amplification size\n ");
	fprintf( bmdata,"%14.6e\n",amplify_step0);

	fclose( bmdata );

  	return 1;    /* ANSI C requires main to return int. */
}

