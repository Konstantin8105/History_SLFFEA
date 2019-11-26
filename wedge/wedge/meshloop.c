/*
    meshloop.c
    This program generates the mesh for a loop based on wedge elements.

     Updated 12/9/08

    SLFFEA source file
    Version:  1.5
    Copyright (C) 2008  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
 */

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi              3.141592654

int mesh_loop(void)
{
	FILE *o1, *o2;
	int i, j, j2, dum;
	int numnp, numel, numel_cross_sec, numel_loop_arc, el_cross_sec, numnp_cross_sec;
	double angle, S, C, Ro, arc_R, RT, X, Y, Z;
	double arc_angle, c_arc_angle, cos_arc_ang, sin_arc_ang,
		c_loop_angle, dloop_angle, darc_angle;
	o1 = fopen( "loop", "w" );
	o2 = fopen( "loop.in", "r" );
	if(o2 == NULL ) {
		printf("Can't find file loop.in\n");
		exit(1);
	}

	char name[30], buf[ BUFSIZ ];

	fgets( buf, BUFSIZ, o2 );
	fscanf( o2, "%d\n",&numel_cross_sec);  /* elements in the loop's cross section */
	printf( "elements in the loop's cross section: %d\n",numel_cross_sec);
	fgets( buf, BUFSIZ, o2 );
	fscanf( o2, "%lf\n",&Ro);              /* loop Radius */
	printf( "loop outer Radius: %lf\n",Ro);
	fgets( buf, BUFSIZ, o2 );
	fscanf( o2, "%lf\n",&arc_R);           /* loop arc Radius */
	printf( "loop arc Radius: %lf\n",arc_R);
	fgets( buf, BUFSIZ, o2 );
	fscanf( o2, "%lf\n",&arc_angle);       /* loop arc Angle */
	printf( "loop arc Angle: %lf\n",arc_angle);
	fgets( buf, BUFSIZ, o2 );
	fscanf( o2, "%d\n",&numel_loop_arc);  /* elements in the loop arc */
	printf( "elements in the loop arc: %d\n",numel_loop_arc);

	numnp_cross_sec = numel_cross_sec + 1;
	numnp = numnp_cross_sec*(numel_loop_arc+1);
	numel = numel_cross_sec*numel_loop_arc;

	fprintf(o1, "   numel numnp nmat nmode   (This is for a wedge loop)\n");
	fprintf(o1, "%d %d %d %d\n", numel, numnp, 1, 0);
	fprintf(o1, "matl no., E modulus, Poisson Ratio, density\n");
	fprintf(o1, "%d %16.4f %16.4f %16.4f\n",0, 1000000.0, 0.27, 2.77e3);
	fprintf(o1, "el no.,connectivity, matl no.\n");
	dum = 0;
	for( i = 0; i < numel_loop_arc; ++i )
	{
	    for( j = 0; j < numel_cross_sec-1; ++j )
	    {
		fprintf(o1, "%4d   ",dum);
		++dum;
		fprintf(o1, "%4d %4d %4d %4d %4d %4d %4d\n",
			numnp_cross_sec*i,
			j+2 + numnp_cross_sec*i, j+1 + numnp_cross_sec*i,
			numnp_cross_sec*(i+1),
			j+2 + numnp_cross_sec*(i+1), j+1 + numnp_cross_sec*(i+1),
			0);
	    }
	    fprintf(o1, "%4d   ",dum);
	    ++dum;
	    j = numel_cross_sec-1;
	    fprintf(o1, "%4d %4d %4d %4d %4d %4d %4d\n",
		numnp_cross_sec*i,
		1 + numnp_cross_sec*i, j+1 + numnp_cross_sec*i,
		numnp_cross_sec*(i+1),
		1 + numnp_cross_sec*(i+1), j+1 + numnp_cross_sec*(i+1),
		0);
	}
	fprintf(o1, "node no. coordinates\n");

	c_arc_angle=0.0;
	dum = 0;
	dloop_angle = 2*pi/(double)numel_cross_sec;
	darc_angle = (arc_angle*pi/180.0)/(double)numel_loop_arc;

	c_arc_angle = 0.0;
	for( i = 0; i < numel_loop_arc + 1; ++i )
	{
	    cos_arc_ang = cos(c_arc_angle);
	    sin_arc_ang = sin(c_arc_angle);
	    c_loop_angle=0.0;

	    X = arc_R*sin_arc_ang;
	    Y = 0.0;
	    Z = arc_R*cos_arc_ang;
	    fprintf(o1, "%4d ",dum);
	    ++dum;
	    fprintf(o1, "%9.5f %9.5f %9.5f \n", X, Y, Z);

	    for( j = 0; j < numel_cross_sec; ++j )
	    {
		X = Ro*cos(c_loop_angle)*sin_arc_ang + arc_R*sin_arc_ang;
		Y = Ro*sin(c_loop_angle);
		Z = Ro*cos(c_loop_angle)*cos_arc_ang + arc_R*cos_arc_ang;
		fprintf(o1, "%4d ",dum);
		++dum;
		fprintf(o1, "%9.5f %9.5f %9.5f \n", X, Y, Z);

		c_loop_angle += dloop_angle;
	    }
	    c_arc_angle += darc_angle;
	}
	fprintf( o1, "prescribed displacement x: node  disp value\n");
	for( j = 0; j < numnp_cross_sec; ++j )
	{
		fprintf(o1, "%4d %9.5f \n", j, 0.0);
	}
	fprintf( o1, "%4d\n ",-10);
	fprintf( o1, "prescribed displacement y: node  disp value\n");
	for( j = 0; j < numnp_cross_sec; ++j )
	{
		fprintf(o1, "%4d %9.5f \n", j, 0.0);
	}
	fprintf( o1, "%4d\n ",-10);
	fprintf( o1, "prescribed displacement z: node  disp value\n");
	for( j = 0; j < numnp_cross_sec; ++j )
	{
		fprintf(o1, "%4d %9.5f \n", j, 0.0);
	}
	fprintf( o1, "%4d\n ",-10);
	fprintf( o1, "node with point load and load vector in x,y,z\n");
	for( j = 0; j < numnp_cross_sec; ++j )
	{
		fprintf(o1, "%4d %9.5f  %9.5f  %9.5f\n",
			j+numel_loop_arc*numnp_cross_sec, -1000.0, 0.0, 0.0);
	}
	fprintf( o1, "%4d\n ",-10);
	fprintf( o1, "node no. with stress and stress vector xx,yy,xy,zx,yz\n");
	fprintf( o1, "%4d ",-10);

	return 1;
}

main(int argc, char** argv)
{
	int check;
	check=mesh_loop();
} 
