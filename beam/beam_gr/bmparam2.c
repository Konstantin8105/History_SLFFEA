/*
    This program calculates and writes the parameters for
    the FEM GUI for beam elements.
  
   			Last Update 4/22/00

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
#include "../beam/bmconst.h"
#include "../beam/bmstruct.h"
#include "bmstrcgr.h"
#include "../../common_gr/control.h"

#define init_far0      -2.0

extern int nmat, numnp, numel, dof;
extern double step_sizex, step_sizey, step_sizez;
extern double left, right, top, bottom, near, far, fscale;
extern int control_height, control_width, mesh_height, mesh_width;
extern double ortho_left, ortho_right, ortho_top, ortho_bottom;
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

void bmReGetparameter(void)
{
        int i, j, check; 
	int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min, node_Uz_max, node_Uz_min,
		node_Uphi_x_max, node_Uphi_x_min, node_Uphi_y_max, node_Uphi_y_min,
		node_Uphi_z_max, node_Uphi_z_min;
	IMDIM max_moment_el, min_moment_el, max_moment_integ, min_moment_integ,
		max_curve_el, min_curve_el, max_curve_integ, min_curve_integ;
	ISDIM max_stress_el, min_stress_el, max_stress_integ, min_stress_integ,
		max_strain_el, min_strain_el, max_strain_integ, min_strain_integ;
	char char_dum[20], char_dum2[5], char_dum3[5], buf[ BUFSIZ ];
	double fdum;
	FILE *bmdata;

/*   bmdata contains all the parameters and extreme values
*/
        bmdata = fopen( "bmview.dat","r" );

/* Read data from the file "bmview.dat" */

	fgets( buf, BUFSIZ, bmdata );
	fgets( buf, BUFSIZ, bmdata );
	fscanf( bmdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Ux_min,
		&node_Ux_max, &min_Ux, &max_Ux);
	fscanf( bmdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Uy_min,
		&node_Uy_max, &min_Uy, &max_Uy);
	fscanf( bmdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Uz_min,
		&node_Uz_max, &min_Uz, &max_Uz);
	fscanf( bmdata,"%20s %5s %5s  %d %d   %lf %lf\n", char_dum, char_dum2, char_dum3,
		&node_Uphi_x_min, &node_Uphi_x_max, &min_Uphi_x, &max_Uphi_x);
	fscanf( bmdata,"%20s %5s %5s  %d %d   %lf %lf\n", char_dum, char_dum2, char_dum3,
		&node_Uphi_y_min, &node_Uphi_y_max, &min_Uphi_y, &max_Uphi_y);
	fscanf( bmdata,"%20s %5s %5s  %d %d   %lf %lf\n", char_dum, char_dum2, char_dum3,
		&node_Uphi_z_min, &node_Uphi_z_max, &min_Uphi_z, &max_Uphi_z);
	fscanf( bmdata,"\n");
	fgets( buf, BUFSIZ, bmdata );
	fgets( buf, BUFSIZ, bmdata );
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_el.xx, &min_moment_integ.xx, &max_moment_el.xx,
		&max_moment_integ.xx, &min_moment.xx, &max_moment.xx);
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_el.yy, &min_moment_integ.yy, &max_moment_el.yy,
		&max_moment_integ.yy, &min_moment.yy, &max_moment.yy);
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_el.zz, &min_moment_integ.zz, &max_moment_el.zz,
		&max_moment_integ.zz, &min_moment.zz, &max_moment.zz);
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_el.xx, &min_stress_integ.xx, &max_stress_el.xx,
		&max_stress_integ.xx, &min_stress.xx, &max_stress.xx);
	fscanf( bmdata,"\n");
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_el.xx, &min_curve_integ.xx, &max_curve_el.xx,
		&max_curve_integ.xx, &min_curve.xx, &max_curve.xx);
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_el.yy, &min_curve_integ.yy, &max_curve_el.yy,
		&max_curve_integ.yy, &min_curve.yy, &max_curve.yy);
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_el.zz, &min_curve_integ.zz, &max_curve_el.zz,
		&max_curve_integ.zz, &min_curve.zz, &max_curve.zz);
	fscanf( bmdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_el.xx, &min_strain_integ.xx, &max_strain_el.xx,
		&max_strain_integ.xx, &min_strain.xx, &max_strain.xx);
	fscanf( bmdata,"\n");
	fgets( buf, BUFSIZ, bmdata );
	fscanf( bmdata,"%lf %lf %lf %lf %lf %lf\n", &ortho_right, &ortho_left,
		&ortho_top, &ortho_bottom, &near, &fdum);
	fgets( buf, BUFSIZ, bmdata );
	fscanf( bmdata,"%6d %6d\n", &mesh_width, &mesh_height);
	fgets( buf, BUFSIZ, bmdata );
	fscanf( bmdata,"%lf %lf %lf\n", &step_sizex, &step_sizey, &step_sizez);
	fgets( buf, BUFSIZ, bmdata );
	fscanf( bmdata,"%lf\n", &amplify_step0);

	fclose( bmdata );

	printf( "                            node\n");
	printf( "                          min  max       min            max\n");
	printf("displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux, max_Ux);
	printf("displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy, max_Uy);
	printf("displacement Uz        %5d %5d   %14.6e %14.6e\n", node_Uz_min,
		node_Uz_max, min_Uz, max_Uz);
	printf("angle phi x            %5d %5d   %14.6e %14.6e\n", node_Uphi_x_min,
		node_Uphi_x_max, min_Uphi_x, max_Uphi_x);
	printf("angle phi y            %5d %5d   %14.6e %14.6e\n", node_Uphi_y_min,
		node_Uphi_y_max, min_Uphi_y, max_Uphi_y);
	printf("angle phi z            %5d %5d   %14.6e %14.6e\n", node_Uphi_z_min,
		node_Uphi_z_max, min_Uphi_z, max_Uphi_z);
	printf("\n");
	printf( "                        el. gauss pt.\n");
	printf( "                        min       max         min           max\n");
	printf("moment xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_moment_el.xx,
		min_moment_integ.xx, max_moment_el.xx, max_moment_integ.xx,
		min_moment.xx, max_moment.xx);
	printf("moment yy            %5d %2d %5d %2d  %14.6e %14.6e\n", min_moment_el.yy,
		min_moment_integ.yy, max_moment_el.yy, max_moment_integ.yy,
		min_moment.yy, max_moment.yy);
	printf("moment zz            %5d %2d %5d %2d  %14.6e %14.6e\n", min_moment_el.zz,
		min_moment_integ.zz, max_moment_el.zz, max_moment_integ.zz,
		min_moment.zz, max_moment.zz);
	printf("stress xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_stress_el.xx,
		min_stress_integ.xx, max_stress_el.xx, max_stress_integ.xx,
		min_stress.xx, max_stress.xx);
	printf("\n");
	printf("curve xx             %5d %2d %5d %2d  %14.6e %14.6e\n", min_curve_el.xx,
		min_curve_integ.xx, max_curve_el.xx, max_curve_integ.xx,
		min_curve.xx, max_curve.xx);
	printf("curve yy             %5d %2d %5d %2d  %14.6e %14.6e\n", min_curve_el.yy,
		min_curve_integ.yy, max_curve_el.yy, max_curve_integ.yy,
		min_curve.yy, max_curve.yy);
	printf("curve zz             %5d %2d %5d %2d  %14.6e %14.6e\n", min_curve_el.zz,
		min_curve_integ.zz, max_curve_el.zz, max_curve_integ.zz,
		min_curve.zz, max_curve.zz);
	printf("strain xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_strain_el.xx,
		min_strain_integ.xx, max_strain_el.xx, max_strain_integ.xx,
		min_strain.xx, max_strain.xx);
	printf("\n");
	printf("Orthographic viewport parameters(right, left, top, bootom, near, far)\n ");
	printf("%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ortho_right, ortho_left,
		ortho_top, ortho_bottom, near, 1000.0);
	printf("Perspective viewport parameters( mesh width and height)\n ");
	printf("%6d %6d\n", mesh_width, mesh_height);
	printf("Step sizes in x, y, z\n ");
	printf("%14.6e %14.6e %14.6e\n",step_sizex, step_sizey, step_sizez);
	printf("Amplification size\n ");
	printf("%14.6e\n",amplify_step0);
}

