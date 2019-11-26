/*
    This program calculates and writes the parameters for
    the FEM GUI for plate elements.
  
   			Last Update 1/23/01

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
#include "../plate/plconst.h"
#include "../plate/plstruct.h"
#include "plstrcgr.h"
#include "../../common_gr/control.h"

#define init_far0      -2.0

extern int nmat, numnp, numel, dof;
extern double step_sizex, step_sizey, step_sizez;
extern double left, right, top, bottom, near, far, fscale, coord_rescale;
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
       	max_Uz, min_Uz, del_Uz, absolute_max_U, absolute_max_coord;

void plReGetparameter(void)
{
        int i, j, check;
        int node_Uz_max, node_Uz_min, node_Uphi_x_max, node_Uphi_x_min,
		node_Uphi_y_max, node_Uphi_y_min;
        IMDIM max_moment_node, min_moment_node, max_curve_node, min_curve_node;
        ISDIM max_stress_node, min_stress_node, max_strain_node, min_strain_node;
	char char_dum[20], char_dum2[5], char_dum3[5], buf[ BUFSIZ ];
	double fdum;
	FILE *pldata;

/*   pldata contains all the parameters and extreme values
*/
	pldata = fopen( "plview.dat","r" );



/* Read data from the file "plview.dat" */

	fgets( buf, BUFSIZ, pldata );
	fgets( buf, BUFSIZ, pldata );
	fscanf( pldata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Uz_min,
		&node_Uz_max, &min_Uz, &max_Uz);

/* Rescale the displacement data */

	min_Uz /= coord_rescale;
	max_Uz /= coord_rescale;

	fscanf( pldata,"%20s %5s %5s  %d %d   %lf %lf\n", char_dum, char_dum2, char_dum3,
		&node_Uphi_x_min, &node_Uphi_x_max, &min_Uphi_x, &max_Uphi_x);
	fscanf( pldata,"%20s %5s %5s  %d %d   %lf %lf\n", char_dum, char_dum2, char_dum3,
		&node_Uphi_y_min, &node_Uphi_y_max, &min_Uphi_y, &max_Uphi_y);
	fscanf( pldata,"\n");
	fgets( buf, BUFSIZ, pldata );
	fgets( buf, BUFSIZ, pldata );
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_node.xx, &max_moment_node.xx,
		&min_moment.xx, &max_moment.xx);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_node.yy, &max_moment_node.yy,
		&min_moment.yy, &max_moment.yy);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_node.xy, &max_moment_node.xy,
		&min_moment.xy, &max_moment.xy);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.zx, &max_stress_node.zx,
		&min_stress.zx, &max_stress.zx);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.yz, &max_stress_node.yz,
		&min_stress.yz, &max_stress.yz);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_node.I, &max_moment_node.I,
		&min_moment.I, &max_moment.I);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_moment_node.II, &max_moment_node.II,
		&min_moment.II, &max_moment.II);
	fscanf( pldata,"\n");
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_node.xx, &max_curve_node.xx,
		&min_curve.xx, &max_curve.xx);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_node.yy, &max_curve_node.yy,
		&min_curve.yy, &max_curve.yy);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_node.xy, &max_curve_node.xy,
		&min_curve.xy, &max_curve.xy);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.zx, &max_strain_node.zx,
		&min_strain.zx, &max_strain.zx);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.yz, &max_strain_node.yz,
		&min_strain.yz, &max_strain.yz);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_node.I, &max_curve_node.I,
		&min_curve.I, &max_curve.I);
	fscanf( pldata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_curve_node.II, &max_curve_node.II,
		&min_curve.II, &max_curve.II);
	fscanf( pldata,"\n");
	fgets( buf, BUFSIZ, pldata );
	fscanf( pldata,"%lf %lf %lf %lf %lf %lf\n", &ortho_right, &ortho_left,
		&ortho_top, &ortho_bottom, &near, &fdum);
	fgets( buf, BUFSIZ, pldata );
	fscanf( pldata,"%d %d\n", &mesh_width, &mesh_height);
	fgets( buf, BUFSIZ, pldata );
	fscanf( pldata,"%lf %lf %lf\n", &step_sizex, &step_sizey, &step_sizez);
	fgets( buf, BUFSIZ, pldata );
	fscanf( pldata,"%lf\n", &amplify_step0);

	fclose( pldata );

	printf( "                            node\n");
	printf( "                          min  max       min            max\n");
	printf( "displacement Uz        %5d %5d   %14.6e %14.6e\n", node_Uz_min,
		node_Uz_max, min_Uz*coord_rescale, max_Uz*coord_rescale);
	printf( "angle phi x            %5d %5d   %14.6e %14.6e\n", node_Uphi_x_min,
		node_Uphi_x_max, min_Uphi_x, max_Uphi_x);
	printf( "angle phi y            %5d %5d   %14.6e %14.6e\n", node_Uphi_y_min,
		node_Uphi_y_max, min_Uphi_y, max_Uphi_y);
	printf( "\n");
	printf( "                            node\n");
	printf( "                        min       max         min           max\n");
	printf( "moment xx            %5d     %5d   %14.6e %14.6e\n", min_moment_node.xx,
		max_moment_node.xx, min_moment.xx, max_moment.xx);
	printf( "moment yy            %5d     %5d   %14.6e %14.6e\n", min_moment_node.yy,
		max_moment_node.yy, min_moment.yy, max_moment.yy);
	printf( "moment xy            %5d     %5d   %14.6e %14.6e\n", min_moment_node.xy,
		max_moment_node.xy, min_moment.xy, max_moment.xy);
	printf( "stress zx            %5d     %5d   %14.6e %14.6e\n", min_stress_node.zx,
		max_stress_node.zx, min_stress.zx, max_stress.zx);
	printf( "stress yz            %5d     %5d   %14.6e %14.6e\n", min_stress_node.yz,
		max_stress_node.yz, min_stress.yz, max_stress.yz);
	printf( "moment I             %5d     %5d   %14.6e %14.6e\n", min_moment_node.I,
		max_moment_node.I, min_moment.I, max_moment.I);
	printf( "moment II            %5d     %5d   %14.6e %14.6e\n", min_moment_node.II,
		max_moment_node.II, min_moment.II, max_moment.II);
	printf( "\n");
	printf( "curve xx             %5d     %5d   %14.6e %14.6e\n", min_curve_node.xx,
		max_curve_node.xx, min_curve.xx, max_curve.xx);
	printf( "curve yy             %5d     %5d   %14.6e %14.6e\n", min_curve_node.yy,
		max_curve_node.yy, min_curve.yy, max_curve.yy);
	printf( "curve xy             %5d     %5d   %14.6e %14.6e\n", min_curve_node.xy,
		max_curve_node.xy, min_curve.xy, max_curve.xy);
	printf( "strain zx            %5d     %5d   %14.6e %14.6e\n", min_strain_node.zx,
		max_strain_node.zx, min_strain.zx, max_strain.zx);
	printf( "strain yz            %5d     %5d   %14.6e %14.6e\n", min_strain_node.yz,
		max_strain_node.yz, min_strain.yz, max_strain.yz);
	printf( "curve I              %5d     %5d   %14.6e %14.6e\n", min_curve_node.I,
		max_curve_node.I, min_curve.I, max_curve.I);
	printf( "curve II             %5d     %5d   %14.6e %14.6e\n", min_curve_node.II,
		max_curve_node.II, min_curve.II, max_curve.II);
	printf( "\n");
	printf( "Orthographic viewport parameters(right, left, top, bootom, near, far)\n ");
	printf( "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n", ortho_right, ortho_left,
		ortho_top, ortho_bottom, near, 1000.0);
	printf( "Perspective viewport parameters( mesh width and height)\n ");
	printf( "%6d %6d\n", mesh_width, mesh_height);
	printf( "Step sizes in x, y, z\n ");
	printf( "%14.6e %14.6e %14.6e\n",step_sizex, step_sizey, step_sizez);
	printf( "Amplification size\n ");
	printf( "%14.6e\n",amplify_step0);

}

