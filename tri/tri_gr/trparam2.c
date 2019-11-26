/*
    This program calculates and writes the parameters for
    the FEM GUI for triangle elements.
  
   			Last Update 6/4/00

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
#include "../tri/trconst.h"
#include "../tri/trstruct.h"
#include "trstrcgr.h"
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
extern STRESS del_stress, del_strain, max_stress, min_stress,
		max_strain, min_strain;
extern double max_Ux, min_Ux, del_Ux, max_Uy, min_Uy, del_Uy,
	       	max_Uz, min_Uz, del_Uz, absolute_max_U;

void trReGetparameter(void)
{
        int i, j, check;
	int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min;
	ISTRESS max_stress_node, min_stress_node, max_strain_node, min_strain_node;
	char char_dum[20], char_dum2[5], buf[ BUFSIZ ];
	double fdum;
	FILE *trdata;
	
/*   trdata contains all the parameters and extreme values
*/
	trdata = fopen( "trview.dat","r" );

	fgets( buf, BUFSIZ, trdata );
	fgets( buf, BUFSIZ, trdata );
	fscanf( trdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Ux_min,
		&node_Ux_max, &min_Ux, &max_Ux);
	fscanf( trdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Uy_min,
		&node_Uy_max, &min_Uy, &max_Uy);

/* Rescale the displacement data */

	min_Ux /= coord_rescale;
	max_Ux /= coord_rescale;
	min_Uy /= coord_rescale;
	max_Uy /= coord_rescale;

	fscanf( trdata,"\n");
	fgets( buf, BUFSIZ, trdata );
	fgets( buf, BUFSIZ, trdata );
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.xx, &max_stress_node.xx,
		&min_stress.xx, &max_stress.xx);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.yy, &max_stress_node.yy,
		&min_stress.yy, &max_stress.yy);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.xy, &max_stress_node.xy,
		&min_stress.xy, &max_stress.xy);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.I, &max_stress_node.I,
		&min_stress.I, &max_stress.I);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.II, &max_stress_node.II,
		&min_stress.II, &max_stress.II);
	fscanf( trdata,"\n");
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.xx, &max_strain_node.xx,
		&min_strain.xx, &max_strain.xx);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.yy, &max_strain_node.yy,
		&min_strain.yy, &max_strain.yy);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.xy, &max_strain_node.xy,
		&min_strain.xy, &max_strain.xy);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.I, &max_strain_node.I,
		&min_strain.I, &max_strain.I);
	fscanf( trdata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.II, &max_strain_node.II,
		&min_strain.II, &max_strain.II);
	fscanf( trdata,"\n");
	fgets( buf, BUFSIZ, trdata );
	fscanf( trdata,"%lf %lf %lf %lf %lf %lf\n", &ortho_right, &ortho_left,
		&ortho_top, &ortho_bottom, &near, &fdum);
	fgets( buf, BUFSIZ, trdata );
	fscanf( trdata,"%d %d\n", &mesh_width, &mesh_height);
	fgets( buf, BUFSIZ, trdata );
	fscanf( trdata,"%lf %lf %lf\n", &step_sizex, &step_sizey, &step_sizez);
	fgets( buf, BUFSIZ, trdata );
	fscanf( trdata,"%lf\n", &amplify_step0);

	fclose( trdata );

	printf( "                            node\n");
	printf( "                          min  max       min            max\n");
	printf("displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux*coord_rescale, max_Ux*coord_rescale);
	printf("displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy*coord_rescale, max_Uy*coord_rescale);
	printf("\n");
	printf( "                        el. gauss pt.\n");
	printf( "                        min       max         min           max\n");
	printf("stress xx            %5d     %5d  %14.6e %14.6e\n", min_stress_node.xx,
		max_stress_node.xx, min_stress.xx, max_stress.xx);
	printf("stress yy            %5d     %5d  %14.6e %14.6e\n", min_stress_node.yy,
		max_stress_node.yy, min_stress.yy, max_stress.yy);
	printf("stress xy            %5d     %5d  %14.6e %14.6e\n", min_stress_node.xy,
		max_stress_node.xy, min_stress.xy, max_stress.xy);
	printf("stress I             %5d     %5d  %14.6e %14.6e\n", min_stress_node.I,
		max_stress_node.I, min_stress.I, max_stress.I);
	printf("stress II            %5d     %5d  %14.6e %14.6e\n", min_stress_node.II,
		max_stress_node.II, min_stress.II, max_stress.II);
	printf("\n");
	printf("strain xx            %5d     %5d  %14.6e %14.6e\n", min_strain_node.xx,
		max_strain_node.xx, min_strain.xx, max_strain.xx);
	printf("strain yy            %5d     %5d  %14.6e %14.6e\n", min_strain_node.yy,
		max_strain_node.yy, min_strain.yy, max_strain.yy);
	printf("strain xy            %5d     %5d  %14.6e %14.6e\n", min_strain_node.xy,
		max_strain_node.xy, min_strain.xy, max_strain.xy);
	printf("strain I             %5d     %5d  %14.6e %14.6e\n", min_strain_node.I,
		max_strain_node.I, min_strain.I, max_strain.I);
	printf("strain II            %5d     %5d  %14.6e %14.6e\n", min_strain_node.II,
		max_strain_node.II, min_strain.II, max_strain.II);
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

