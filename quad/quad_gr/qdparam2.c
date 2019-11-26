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
extern SDIM del_stress, del_strain, max_stress, min_stress,
	max_strain, min_strain;
extern double max_Ux, min_Ux, del_Ux, max_Uy, min_Uy, del_Uy,
	max_Uz, min_Uz, del_Uz, absolute_max_U;

void qdReGetparameter(void)
{
        int i, j, check;
	int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min;
	ISDIM max_stress_node, min_stress_node, max_strain_node, min_strain_node;
	char char_dum[20], char_dum2[5], buf[ BUFSIZ ];
	double fdum;
	FILE *qddata;
	
/*   qddata contains all the parameters and extreme values
*/
	qddata = fopen( "qdview.dat","r" );

	fgets( buf, BUFSIZ, qddata );
	fgets( buf, BUFSIZ, qddata );
	fscanf( qddata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Ux_min,
		&node_Ux_max, &min_Ux, &max_Ux);
	fscanf( qddata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Uy_min,
		&node_Uy_max, &min_Uy, &max_Uy);
	fscanf( qddata,"\n");
	fgets( buf, BUFSIZ, qddata );
	fgets( buf, BUFSIZ, qddata );
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.xx, &max_stress_node.xx,
		&min_stress.xx, &max_stress.xx);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.yy, &max_stress_node.yy,
		&min_stress.yy, &max_stress.yy);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.xy, &max_stress_node.xy,
		&min_stress.xy, &max_stress.xy);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.I, &max_stress_node.I,
		&min_stress.I, &max_stress.I);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_node.II, &max_stress_node.II,
		&min_stress.II, &max_stress.II);
	fscanf( qddata,"\n");
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.xx, &max_strain_node.xx,
		&min_strain.xx, &max_strain.xx);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.yy, &max_strain_node.yy,
		&min_strain.yy, &max_strain.yy);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.xy, &max_strain_node.xy,
		&min_strain.xy, &max_strain.xy);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.I, &max_strain_node.I,
		&min_strain.I, &max_strain.I);
	fscanf( qddata,"%20s %5s    %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_node.II, &max_strain_node.II,
		&min_strain.II, &max_strain.II);
	fscanf( qddata,"\n");
	fgets( buf, BUFSIZ, qddata );
	fscanf( qddata,"%lf %lf %lf %lf %lf %lf\n", &ortho_right, &ortho_left,
		&ortho_top, &ortho_bottom, &near, &fdum);
	fgets( buf, BUFSIZ, qddata );
	fscanf( qddata,"%d %d\n", &mesh_width, &mesh_height);
	fgets( buf, BUFSIZ, qddata );
	fscanf( qddata,"%lf %lf %lf\n", &step_sizex, &step_sizey, &step_sizez);
	fgets( buf, BUFSIZ, qddata );
	fscanf( qddata,"%lf\n", &amplify_step0);

	fclose( qddata );

	printf( "                            node\n");
	printf( "                          min  max       min            max\n");
	printf("displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux, max_Ux);
	printf("displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy, max_Uy);
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

