/*
    This program calculates and writes the parameters for
    the FEM GUI for truss elements.
  
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
#include "../truss/tsconst.h"
#include "../truss/tsstruct.h"
#include "tsstrcgr.h"
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
extern STRESS del_stress, max_stress, min_stress;
extern STRAIN del_strain, max_strain, min_strain;
extern double max_Ux, min_Ux, del_Ux, max_Uy, min_Uy, del_Uy,
	max_Uz, min_Uz, del_Uz, absolute_max_U;

void tsReGetparameter(void)
{
        int i, j, check;
        int node_Ux_max, node_Ux_min, node_Uy_max, node_Uy_min, node_Uz_max, node_Uz_min;
        ISDIM max_stress_el, min_stress_el, max_stress_integ, min_stress_integ,
                max_strain_el, min_strain_el, max_strain_integ, min_strain_integ;
	char char_dum[20], char_dum2[5], buf[ BUFSIZ ];
	double fdum;
        FILE *tsdata;

/*   shdata contains all the parameters and extreme values
*/
        tsdata = fopen( "tsview.dat","r" );

/* Read data from the file "tsview.dat" */

	fgets( buf, BUFSIZ, tsdata );
	fgets( buf, BUFSIZ, tsdata );
	fscanf( tsdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Ux_min,
		&node_Ux_max, &min_Ux, &max_Ux);
	fscanf( tsdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Uy_min,
		&node_Uy_max, &min_Uy, &max_Uy);
	fscanf( tsdata,"%20s %5s      %d %d   %lf %lf\n", char_dum, char_dum2, &node_Uz_min,
		&node_Uz_max, &min_Uz, &max_Uz);
	fscanf( tsdata,"\n");
	fgets( buf, BUFSIZ, tsdata );
	fgets( buf, BUFSIZ, tsdata );
	fscanf( tsdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_stress_el.xx, &min_stress_integ.xx, &max_stress_el.xx,
		&max_stress_integ.xx, &min_stress.xx, &max_stress.xx);
	fscanf( tsdata,"\n");
	fscanf( tsdata,"%20s %5s    %d %d %d %d  %lf %lf\n", char_dum, char_dum2,
		&min_strain_el.xx, &min_strain_integ.xx, &max_strain_el.xx,
		&max_strain_integ.xx, &min_strain.xx, &max_strain.xx);
	fscanf( tsdata,"\n");
	fgets( buf, BUFSIZ, tsdata );
	fscanf( tsdata,"%lf %lf %lf %lf %lf %lf\n", &ortho_right, &ortho_left,
		&ortho_top, &ortho_bottom, &near, &fdum);
	fgets( buf, BUFSIZ, tsdata );
	fscanf( tsdata,"%d %d\n", &mesh_width, &mesh_height);
	fgets( buf, BUFSIZ, tsdata );
	fscanf( tsdata,"%lf %lf %lf\n", &step_sizex, &step_sizey, &step_sizez);
	fgets( buf, BUFSIZ, tsdata );
	fscanf( tsdata,"%lf\n", &amplify_step0);

	fclose( tsdata );

	printf( "                            node\n");
	printf( "                          min  max       min            max\n");
	printf("displacement Ux        %5d %5d   %14.6e %14.6e\n", node_Ux_min,
		node_Ux_max, min_Ux, max_Ux);
	printf("displacement Uy        %5d %5d   %14.6e %14.6e\n", node_Uy_min,
		node_Uy_max, min_Uy, max_Uy);
	printf("displacement Uz        %5d %5d   %14.6e %14.6e\n", node_Uz_min,
		node_Uz_max, min_Uz, max_Uz);
	printf("\n");
	printf( "                        el. gauss pt.\n");
	printf( "                        min       max         min           max\n");
	printf("stress xx            %5d %2d %5d %2d  %14.6e %14.6e\n", min_stress_el.xx,
		min_stress_integ.xx, max_stress_el.xx, max_stress_integ.xx,
		min_stress.xx, max_stress.xx);
	printf("\n");
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

