/*
    This program generates the mesh for a tensile test.  It replaces
    the older version of meshtensile.c which I have renamed meshtensile.c.old.
    The problems with "meshtensile.c.old" is discussed in the file itself.

    My original intention was to simply fix the problems discussed there,
    but I decided to treat this code as an experiment in how to insert nodes
    into a mesh and how to change the nodal data and connectivity after
    that happens.  I will discuss this further at the points where I do
    it in the code.

    I will also add that the original mesh before any nodes are inserted
    should be the same as that of the meshes generated by "meshtorsion.c".

		Updated 2/15/10

    SLFFEA source file
    Version:  1.6
    Copyright (C) 1999-2010  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define nsd             3
#define npel            8
#define pi              3.141592654
#define sq2             1.414213462
#define SMALL           1.e-20

typedef struct {
	int node;
	int index;
} NEWNODE;

int compare_newnode ( const void *x,  const void *y)
{
	NEWNODE *a = (NEWNODE *) x;
	NEWNODE *b = (NEWNODE *) y;
	if ( a[0].node < b[0].node ) return -1;
	if ( a[0].node == b[0].node ) return 0;
	if ( a[0].node > b[0].node ) return 1;
}

int main(int argc, char** argv)
{
	int dof, nmat, numel, numnp, nmode, newnode_counter;
	int i, i2, j, k, dum, dum2, dum3, dum3a, dum4, dum5, dum6, dum7, dum8;
	int *connect;
	NEWNODE *newnode;
	double fdum, fdum1, fdum2, fdum3, fdum4, ux[9], uy[9];
	double angle, height_length, square_length, cylinder_length, ray_angle,
		height_el_length, square_el_length, cylinder_el_length,
		ratio;
	double *coord;
	int height_el_num, half_height_el_num, square_el_num, cylinder_el_num, angle_div,
		angleh_div, nodes_per_layer;
	char name[30];
	char buf[ BUFSIZ ];
	FILE *o1;
	char text;

	o1 = fopen( "ttest","w" );
	printf( "\n What is the length of the inner square?\n ");
	scanf( "%lf",&square_length);
	printf( "\n How many elements on the inner square?\n ");
	scanf( "%d",&square_el_num);
	printf( "\n What is the length of the cylinder height?\n ");
	scanf( "%lf",&height_length);
	printf( "\n How many elements on the cylinder height?\n ");
	scanf( "%d",&height_el_num);
	printf( "\n What is the length of the outer cylinder?\n ");
	scanf( "%lf",&cylinder_length);
	printf( "\n How many elements on the outer cylinder?\n ");
	scanf( "%d",&cylinder_el_num);
	printf( "\n What is the ratio?(ratio > 1)\n ");
	scanf( "%lf",&ratio);

	ray_angle = 45.0/((double)square_el_num);
	angle_div = (int)(90.0/ray_angle);
	angleh_div = (int)(((double)angle_div)/2.0);
	/*nodes_per_layer = (angle_div+1)*(cylinder_el_num+square_el_num+1);*/
	nodes_per_layer = (square_el_num+1)*(square_el_num+1) +
	    (angle_div + 1)*(cylinder_el_num) +
	    2*(square_el_num)*(square_el_num+1) +
	    2*(angle_div)*(cylinder_el_num) +
	    (square_el_num)*(square_el_num) +
	    (angle_div - 1)*(cylinder_el_num);

	nmat = 1; 

	nmode = 0;
	numel = (square_el_num*square_el_num + (angle_div)*(cylinder_el_num))*height_el_num;
	numel *= 4;
	numnp = nodes_per_layer*(height_el_num+1) + (cylinder_el_num-1)*(square_el_num-2);
	coord=(double *)calloc(nsd*numnp,sizeof(double));
	connect=(int *)calloc(numel*npel,sizeof(int));
	dum = (int)(numnp/10);
	newnode=(NEWNODE *)calloc(dum,sizeof(NEWNODE));

	half_height_el_num = (height_el_num - height_el_num%2)/2;
	cylinder_el_length = cylinder_length/((double)cylinder_el_num);
	square_el_length = square_length/((double)square_el_num);
	height_el_length = height_length/((double)height_el_num);
	printf( "    %4d %4d %4d %4d \n ", numel,numnp,nmat,nodes_per_layer);
	printf( "    %10.6f %10.6f %10.6f\n ", square_el_length, height_el_length, cylinder_el_length );

	fprintf( o1, "   numel numnp nmat nmode  This is the Saint Venant mesh \n ");
	fprintf( o1, "    %4d %4d %4d %4d\n ", numel,numnp,nmat,nmode);

	fprintf( o1, "matl no., E modulus, Poisson Ratio, density \n");
	for( i = 0; i < nmat; ++i )
	{
	   fprintf( o1, " %3d",0);
	   fprintf( o1, " %14.6e %14.6e %14.6e\n ", 2.410000e+11, 0.2, 2.77e3);
	}

	ray_angle *= pi/180.0;
	dum = 0;
	dum2 = 0;
	dum3 = (square_el_num+cylinder_el_num+1);
	dum3a = (square_el_num+cylinder_el_num);
	dum4 = 0;
	dum6 = 0;
	dum8 = (square_el_num+1)*(square_el_num+1) +
	    (angle_div + 1)*(cylinder_el_num) +
	    2*(square_el_num)*(square_el_num+1) +
	    2*(angle_div)*(cylinder_el_num);

	for( k = 0; k < height_el_num; ++k )
	{
	    dum5 = 0;

/* This is the 1st quarter */

	    for( i = 0; i < square_el_num; ++i )
	    {
		for( j = 0; j < square_el_num + cylinder_el_num; ++j )
		{
		    *(connect + npel*dum)     = i*dum3 + j + dum2;
		    *(connect + npel*dum + 1) = i*dum3 + j + 1 + dum2;
		    *(connect + npel*dum + 2) = (i+1)*dum3 + j + 1 + dum2;
		    *(connect + npel*dum + 3) = (i+1)*dum3 + j + dum2;
		    *(connect + npel*dum + 4) = i*dum3 + j + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = i*dum3 + j + 1 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = (i+1)*dum3 + j + 1 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = (i+1)*dum3 + j + nodes_per_layer + dum2;
		    ++dum;
		}
	    }
	    dum4 = square_el_num*dum3;
	    for( j = 0; j < square_el_num - 1; ++j )
	    {
		    *(connect + npel*dum)     = dum4 + j + dum2;
		    *(connect + npel*dum + 1) = dum4 + j + 1 + dum2;
		    *(connect + npel*dum + 2) = dum4 + j + dum3 + 1 + dum2;
		    *(connect + npel*dum + 3) = dum4 + j + dum3 + dum2;
		    *(connect + npel*dum + 4) = dum4 + j + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = dum4 + j + 1 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = dum4 + j + dum3 + 1 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = dum4 + j + dum3 + nodes_per_layer + dum2;
		++dum;
	    }
	    dum4 = (1+square_el_num)*dum3;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		for( j = 0; j < square_el_num - 1; ++j )
		{
		    *(connect + npel*dum)     = dum4 + i*square_el_num + j + dum2;
		    *(connect + npel*dum + 1) = dum4 + i*square_el_num + j + 1 + dum2;
		    *(connect + npel*dum + 2) = dum4 + (i+1)*square_el_num + j + 1 + dum2;
		    *(connect + npel*dum + 3) = dum4 + (i+1)*square_el_num + j + dum2;
		    *(connect + npel*dum + 4) = dum4 + i*square_el_num + j + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = dum4 + i*square_el_num + j + 1 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = dum4 + (i+1)*square_el_num + j + 1 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = dum4 + (i+1)*square_el_num + j + nodes_per_layer + dum2;
		    ++dum;
		}
	    }
	    dum4 = square_el_num*dum3;
	    *(connect + npel*dum)     = dum4 + square_el_num - 1 + dum2;
	    *(connect + npel*dum + 1) = dum4 + square_el_num + dum2;
	    *(connect + npel*dum + 2) = dum4 + square_el_num + 1 + dum2;
	    *(connect + npel*dum + 3) = dum4 + dum3 - 1 + square_el_num + dum2;
	    *(connect + npel*dum + 4) = dum4 + square_el_num - 1 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 5) = dum4 + square_el_num + nodes_per_layer + dum2;
	    *(connect + npel*dum + 6) = dum4 + square_el_num + 1 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 7) = dum4 + dum3 - 1 + square_el_num + nodes_per_layer + dum2;
	    ++dum;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		*(connect + npel*dum)    = (i+1)*square_el_num - 1 + dum3 + dum4 + dum2;
		*(connect + npel*dum + 1) = dum4 + square_el_num + i + 1 + dum2;
		*(connect + npel*dum + 2) = dum4 + square_el_num + i + 2 + dum2;
		*(connect + npel*dum + 3) = (i+1)*square_el_num - 1 + dum3 + dum4 + square_el_num + dum2;
		*(connect + npel*dum + 4) = (i+1)*square_el_num - 1 + dum3 + dum4 + nodes_per_layer + dum2;
		*(connect + npel*dum + 5) = dum4 + square_el_num + i + 1 + nodes_per_layer + dum2;
		*(connect + npel*dum + 6) = dum4 + square_el_num + i + 2 + nodes_per_layer + dum2;
		*(connect + npel*dum + 7) =
			(i+1)*square_el_num - 1 + dum3 + dum4 + square_el_num + nodes_per_layer + dum2;
		++dum;
	    }

/* This is the 2nd quarter */

	    dum5 += (square_el_num+1)*(square_el_num+1) + (angle_div+1)*(cylinder_el_num);

	    for( j = 0; j < square_el_num + 1; ++j )
	    {
		    *(connect + npel*dum)     = j*(square_el_num + cylinder_el_num + 1) + dum2;
		    *(connect + npel*dum + 1) = (j + 1)*(square_el_num + cylinder_el_num + 1) + dum2;
		    *(connect + npel*dum + 2) = j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = j + dum5 + dum2;
		    *(connect + npel*dum + 4) = j*(square_el_num + cylinder_el_num + 1) + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = (j + 1)*(square_el_num + cylinder_el_num + 1) + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = j + 1 + dum5 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = j + dum5 + nodes_per_layer + dum2;
		    ++dum;
	    }

	    dum6 = (square_el_num + 1)*(square_el_num + cylinder_el_num + 1);
	    for( j = 1; j < cylinder_el_num; ++j )
	    {
		    *(connect + npel*dum)     = (j - 1)*(square_el_num) + dum6 + dum2;
		    *(connect + npel*dum + 1) = j*square_el_num + dum6 + dum2;
		    *(connect + npel*dum + 2) = j + square_el_num + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = j + square_el_num + dum5 + dum2;
		    *(connect + npel*dum + 4) = (j - 1)*(square_el_num) + dum6 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = j*square_el_num + dum6 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = j + square_el_num + 1 + dum5 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = j + square_el_num + dum5 + nodes_per_layer + dum2;
		    ++dum;
	    }

	    for( i = 1; i < square_el_num; ++i )
	    {
		for( j = 0; j < square_el_num + cylinder_el_num; ++j )
		{
		    *(connect + npel*dum)     = (i-1)*dum3 + j + dum5 + dum2;
		    *(connect + npel*dum + 1) = (i-1)*dum3 + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 2) = i*dum3 + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = i*dum3 + j + dum5 + dum2;
		    *(connect + npel*dum + 4) = (i-1)*dum3 + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 5) = (i-1)*dum3 + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 6) = i*dum3 + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 7) = i*dum3 + j + nodes_per_layer + dum5 + dum2;
		    ++dum;
		}
	    }
	    dum4 = square_el_num*dum3;
	    for( j = 0; j < square_el_num - 1; ++j )
	    {
		*(connect + npel*dum)     = dum3*(square_el_num - 1) + j + dum5 + dum2;
		*(connect + npel*dum + 1) = dum3*(square_el_num - 1) + j + 1 + dum5 + dum2;
		*(connect + npel*dum + 2) = dum3*square_el_num + j + 1 + dum5 + dum2;
		*(connect + npel*dum + 3) = dum3*square_el_num + j + dum5 + dum2;
		*(connect + npel*dum + 4) = dum3*(square_el_num - 1) + j + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 5) = dum3*(square_el_num - 1) + j + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 6) = dum3*square_el_num + j + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 7) = dum3*square_el_num + j + nodes_per_layer + dum5 + dum2;
		++dum;
	    }
	    dum4 = (1+square_el_num)*dum3;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		for( j = 0; j < square_el_num - 1; ++j )
		{
		    *(connect + npel*dum)     = dum3*square_el_num + i*square_el_num + j + dum5 + dum2;
		    *(connect + npel*dum + 1) = dum3*square_el_num + i*square_el_num + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 2) = dum3*square_el_num + (i+1)*square_el_num + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = dum3*square_el_num + (i+1)*square_el_num + j + dum5 + dum2;
		    *(connect + npel*dum + 4) = dum3*square_el_num + i*square_el_num + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 5) = dum3*square_el_num + i*square_el_num + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 6) = dum3*square_el_num + (i+1)*square_el_num + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 7) = dum3*square_el_num + (i+1)*square_el_num + j + nodes_per_layer + dum5 + dum2;
		    ++dum;
		}
	    }
	    dum4 = square_el_num*dum3;
	    *(connect + npel*dum)     = dum3*(square_el_num-1) + square_el_num - 1 + dum5 + dum2;
	    *(connect + npel*dum + 1) = dum3*(square_el_num-1) + square_el_num + dum5 + dum2;
	    *(connect + npel*dum + 2) = dum3*(square_el_num-1) + square_el_num + 1 + dum5 + dum2;
	    *(connect + npel*dum + 3) = dum3*square_el_num - 1 + square_el_num + dum5 + dum2;
	    *(connect + npel*dum + 4) = dum3*(square_el_num-1) + square_el_num - 1 + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 5) = dum3*(square_el_num-1) + square_el_num + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 6) = dum3*(square_el_num-1) + square_el_num + 1 + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 7) = dum3*square_el_num - 1 + square_el_num + nodes_per_layer + dum5 + dum2;
	    ++dum;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		*(connect + npel*dum)     = (i+1+dum3)*square_el_num - 1 + dum5 + dum2;
		*(connect + npel*dum + 1) = dum3*(square_el_num-1) + square_el_num + i + 1 + dum5 + dum2;
		*(connect + npel*dum + 2) = dum3*(square_el_num-1) + square_el_num + i + 2 + dum5 + dum2;
		*(connect + npel*dum + 3) = (i+1+dum3)*square_el_num - 1 + square_el_num + dum5 + dum2;
		*(connect + npel*dum + 4) = (i+1+dum3)*square_el_num - 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 5) = dum3*(square_el_num-1) + square_el_num + i + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 6) = dum3*(square_el_num-1) + square_el_num + i + 2 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 7) = (i+1+dum3)*square_el_num - 1 + square_el_num + nodes_per_layer + dum5 + dum2;
		++dum;
	    }

/* This is the 3rd quarter */

	    /*dum7 = dum5;*/
	    dum7 = (square_el_num)*(square_el_num+1) + (angle_div)*(cylinder_el_num);
	    dum5 += (square_el_num)*(square_el_num+1) + (angle_div)*(cylinder_el_num);

	    *(connect + npel*dum)     = dum2;
	    *(connect + npel*dum + 1) = (square_el_num + cylinder_el_num + 1) + dum7 + dum2;
	    *(connect + npel*dum + 2) = 1 + dum5 + dum2;
	    *(connect + npel*dum + 3) = dum5 + dum2;
	    *(connect + npel*dum + 4) = nodes_per_layer + dum2;
	    *(connect + npel*dum + 5) = (square_el_num + cylinder_el_num + 1) + dum7 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 6) = 1 + dum5 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 7) = dum5 + nodes_per_layer + dum2;
	    ++dum;

	    for( j = 1; j < square_el_num + 1; ++j )
	    {
		    *(connect + npel*dum)     = j*(square_el_num + cylinder_el_num + 1) + dum7 + dum2;
		    *(connect + npel*dum + 1) = (j + 1)*(square_el_num + cylinder_el_num + 1) + dum7 + dum2;
		    *(connect + npel*dum + 2) = j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = j + dum5 + dum2;
		    *(connect + npel*dum + 4) = j*(square_el_num + cylinder_el_num + 1) + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = (j + 1)*(square_el_num + cylinder_el_num + 1) + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = j + 1 + dum5 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = j + dum5 + nodes_per_layer + dum2;
		    ++dum;
	    }

	    dum6 = (square_el_num + 1)*(square_el_num + cylinder_el_num + 1);
	    for( j = 1; j < cylinder_el_num; ++j )
	    {
		    *(connect + npel*dum)     = (j - 1)*(square_el_num) + dum6 + dum7 + dum2;
		    *(connect + npel*dum + 1) = j*square_el_num + dum6 + dum7 + dum2;
		    *(connect + npel*dum + 2) = j + square_el_num + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = j + square_el_num + dum5 + dum2;
		    *(connect + npel*dum + 4) = (j - 1)*(square_el_num) + dum6 + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = j*square_el_num + dum6 + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = j + square_el_num + 1 + dum5 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = j + square_el_num + dum5 + nodes_per_layer + dum2;
		    ++dum;
	    }

	    for( i = 1; i < square_el_num; ++i )
	    {
		for( j = 0; j < square_el_num + cylinder_el_num; ++j )
		{
		    *(connect + npel*dum)     = (i-1)*dum3 + j + dum5 + dum2;
		    *(connect + npel*dum + 1) = (i-1)*dum3 + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 2) = i*dum3 + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = i*dum3 + j + dum5 + dum2;
		    *(connect + npel*dum + 4) = (i-1)*dum3 + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 5) = (i-1)*dum3 + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 6) = i*dum3 + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 7) = i*dum3 + j + nodes_per_layer + dum5 + dum2;
		    ++dum;
		}
	    }
	    dum4 = square_el_num*dum3;
	    for( j = 0; j < square_el_num - 1; ++j )
	    {
		*(connect + npel*dum)     = dum3*(square_el_num-1) + j + dum5 + dum2;
		*(connect + npel*dum + 1) = dum3*(square_el_num-1) + j + 1 + dum5 + dum2;
		*(connect + npel*dum + 2) = dum3*square_el_num + j + 1 + dum5 + dum2;
		*(connect + npel*dum + 3) = dum3*square_el_num + j + dum5 + dum2;
		*(connect + npel*dum + 4) = dum3*(square_el_num-1) + j + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 5) = dum3*(square_el_num-1) + j + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 6) = dum3*square_el_num + j + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 7) = dum3*square_el_num + j + nodes_per_layer + dum5 + dum2;
		++dum;
	    }
	    dum4 = (1+square_el_num)*dum3;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		for( j = 0; j < square_el_num - 1; ++j )
		{
		    *(connect + npel*dum)     = (dum3+i)*square_el_num + j + dum5 + dum2;
		    *(connect + npel*dum + 1) = (dum3+i)*square_el_num + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 2) = (dum3+i+1)*square_el_num + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 3) = (dum3+i+1)*square_el_num + j + dum5 + dum2;
		    *(connect + npel*dum + 4) = (dum3+i)*square_el_num + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 5) = (dum3+i)*square_el_num + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 6) = (dum3+i+1)*square_el_num + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 7) = (dum3+i+1)*square_el_num + j + nodes_per_layer + dum5 + dum2;
		    ++dum;
		}
	    }
	    dum4 = square_el_num*dum3;
	    *(connect + npel*dum)     = dum3*(square_el_num-1) + square_el_num - 1 + dum5 + dum2;
	    *(connect + npel*dum + 1) = dum3*(square_el_num-1) + square_el_num + dum5 + dum2;
	    *(connect + npel*dum + 2) = dum3*(square_el_num-1) + square_el_num + 1 + dum5 + dum2;
	    *(connect + npel*dum + 3) = dum3*square_el_num - 1 + square_el_num + dum5 + dum2;
	    *(connect + npel*dum + 4) = dum3*(square_el_num-1) + square_el_num - 1 + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 5) = dum3*(square_el_num-1) + square_el_num + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 6) = dum3*(square_el_num-1) + square_el_num + 1 + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 7) = dum3*square_el_num - 1 + square_el_num + nodes_per_layer + dum5 + dum2;
	    ++dum;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		*(connect + npel*dum)     = (i+1+dum3)*square_el_num - 1 + dum5 + dum2;
		*(connect + npel*dum + 1) = dum3*(square_el_num-1) + square_el_num + i + 1 + dum5 + dum2;
		*(connect + npel*dum + 2) = dum3*(square_el_num-1) + square_el_num + i + 2 + dum5 + dum2;
		*(connect + npel*dum + 3) = (i+1+dum3)*square_el_num - 1 + square_el_num + dum5 + dum2;
		*(connect + npel*dum + 4) = (i+1+dum3)*square_el_num - 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 5) = dum3*(square_el_num-1) + square_el_num + i + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 6) = dum3*(square_el_num-1) + square_el_num + i + 2 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 7) = (i+1+dum3)*square_el_num - 1 + square_el_num + nodes_per_layer + dum5 + dum2;
		++dum;
	    }


/* This is the 4th quarter */

	    /*dum7 = dum5;
	    dum5 += (square_el_num)*(square_el_num+1) + (angle_div)*(cylinder_el_num);*/
	    dum7 = 2*((square_el_num)*(square_el_num+1) + (angle_div)*(cylinder_el_num));
	    dum5 += (square_el_num)*(square_el_num+1) + (angle_div)*(cylinder_el_num);

	    *(connect + npel*dum)     = dum2;
	    *(connect + npel*dum + 1) = dum3 + dum7 + dum2;
	    *(connect + npel*dum + 2) = dum5 + dum2;
	    *(connect + npel*dum + 3) = 1 + dum2;
	    *(connect + npel*dum + 4) = nodes_per_layer + dum2;
	    *(connect + npel*dum + 5) = dum3 + dum7 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 6) = dum5 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 7) = 1 + nodes_per_layer + dum2;
	    ++dum;

	    for( j = 1; j < square_el_num + 1; ++j )
	    {
		    *(connect + npel*dum)     = j*(square_el_num + cylinder_el_num + 1) + dum7 + dum2;
		    *(connect + npel*dum + 1) = (j + 1)*(square_el_num + cylinder_el_num + 1) + dum7 + dum2;
		    *(connect + npel*dum + 2) = j + dum5 + dum2;
		    *(connect + npel*dum + 3) = j - 1 + dum5 + dum2;
		    *(connect + npel*dum + 4) = j*(square_el_num + cylinder_el_num + 1) + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = (j + 1)*(square_el_num + cylinder_el_num + 1) + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = j + dum5 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = j - 1 + dum5 + nodes_per_layer + dum2;
		    ++dum;
	    }

	    dum6 = (square_el_num + 1)*(square_el_num + cylinder_el_num + 1);
	    for( j = 1; j < cylinder_el_num; ++j )
	    {
		    *(connect + npel*dum)     = (j - 1)*(square_el_num) + dum6 + dum7 + dum2;
		    *(connect + npel*dum + 1) = j*square_el_num + dum6 + dum7 + dum2;
		    *(connect + npel*dum + 2) = j + square_el_num + dum5 + dum2;
		    *(connect + npel*dum + 3) = j - 1 + square_el_num + dum5 + dum2;
		    *(connect + npel*dum + 4) = (j - 1)*(square_el_num) + dum6 + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 5) = j*square_el_num + dum6 + dum7 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 6) = j + square_el_num + dum5 + nodes_per_layer + dum2;
		    *(connect + npel*dum + 7) = j - 1 + square_el_num + dum5 + nodes_per_layer + dum2;
		    ++dum;
	    }

	    for( i = 1; i < square_el_num; ++i )
	    {
		*(connect + npel*dum)     = i + dum2;
		*(connect + npel*dum + 1) = (i-1)*dum3a + dum5 + dum2;
		*(connect + npel*dum + 2) = (i)*dum3a + dum5 + dum2;
		*(connect + npel*dum + 3) = i + 1 + dum2;
		*(connect + npel*dum + 4) = i + nodes_per_layer + dum2;
		*(connect + npel*dum + 5) = (i-1)*dum3a + dum5 + nodes_per_layer + dum2;
		*(connect + npel*dum + 6) = (i)*dum3a + dum5 + nodes_per_layer + dum2;
		*(connect + npel*dum + 7) = i + 1 + nodes_per_layer + dum2;
		++dum;

		for( j = 0; j < square_el_num + cylinder_el_num - 1; ++j )
		{
		    *(connect + npel*dum)     = (i-1)*dum3a + j + dum5 + dum2;
		    *(connect + npel*dum + 1) = (i-1)*dum3a + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 2) = i*dum3a + 1 + j + dum5 + dum2;
		    *(connect + npel*dum + 3) = i*dum3a + j + dum5 + dum2;
		    *(connect + npel*dum + 4) = (i-1)*dum3a + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 5) = (i-1)*dum3a + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 6) = i*dum3a + 1 + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 7) = i*dum3a + j + nodes_per_layer + dum5 + dum2;
		    ++dum;
		}
	    }

	    *(connect + npel*dum)     = square_el_num + dum2;
	    *(connect + npel*dum + 1) = dum3a*(square_el_num-1) + dum5 + dum2;
	    *(connect + npel*dum + 2) = dum3a*square_el_num + dum5 + dum2;
	    *(connect + npel*dum + 3) = square_el_num + 1 + dum2;
	    *(connect + npel*dum + 4) = square_el_num + nodes_per_layer + dum2;
	    *(connect + npel*dum + 5) = dum3a*(square_el_num-1) + dum5 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 6) = dum3a*square_el_num + dum5 + nodes_per_layer + dum2;
	    *(connect + npel*dum + 7) = square_el_num + 1 + nodes_per_layer + dum2;
	    ++dum;

	    dum4 = square_el_num*dum3a;
	    for( j = 1; j < square_el_num - 1; ++j )
	    {
		*(connect + npel*dum)     = dum3a*(square_el_num-1) - 1 + j + dum5 + dum2;
		*(connect + npel*dum + 1) = dum3a*(square_el_num-1) + j + dum5 + dum2;
		*(connect + npel*dum + 2) = dum3a*square_el_num + j + dum5 + dum2;
		*(connect + npel*dum + 3) = dum3a*square_el_num -1 + j + dum5 + dum2;
		*(connect + npel*dum + 4) = dum3a*(square_el_num-1) -1 + j + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 5) = dum3a*(square_el_num-1) + j + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 6) = dum3a*square_el_num + j + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 7) = dum3a*square_el_num - 1 + j + nodes_per_layer + dum5 + dum2;
		++dum;
	    }
	    dum4 = (1+square_el_num)*dum3a;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		*(connect + npel*dum)     = square_el_num + i + 1 + dum2;
		*(connect + npel*dum + 1) = i*(square_el_num-1) + dum3a*(square_el_num) + dum5 + dum2;
		*(connect + npel*dum + 2) = (i+1)*(square_el_num-1) + dum3a*(square_el_num) + dum5 + dum2;
		*(connect + npel*dum + 3) = square_el_num + i + 2 + dum2;
		*(connect + npel*dum + 4) = square_el_num + i + 1 + nodes_per_layer + dum2;
		*(connect + npel*dum + 5) = i*(square_el_num-1) + dum3a*(square_el_num) + dum5 + nodes_per_layer + dum2;
		*(connect + npel*dum + 6) = (i+1)*(square_el_num-1) + dum3a*(square_el_num) + dum5 + nodes_per_layer + dum2;
		*(connect + npel*dum + 7) = square_el_num + i + 2 + nodes_per_layer + dum2;
		++dum;

		for( j = 0; j < square_el_num - 3; ++j )
		{
		    *(connect + npel*dum)     = (dum3a+i)*square_el_num - i + j + dum5 + dum2;
		    *(connect + npel*dum + 1) = (dum3a+i)*square_el_num - i + j + 1 + dum5 + dum2;
		    *(connect + npel*dum + 2) = (dum3a+i+1)*square_el_num - i + j + dum5 + dum2;
		    *(connect + npel*dum + 3) = (dum3a+i+1)*square_el_num - i + j - 1 + dum5 + dum2;
		    *(connect + npel*dum + 4) = (dum3a+i)*square_el_num - i + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 5) = (dum3a+i)*square_el_num - i + j + 1 + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 6) = (dum3a+i+1)*square_el_num - i + j + nodes_per_layer + dum5 + dum2;
		    *(connect + npel*dum + 7) = (dum3a+i+1)*square_el_num - i + j - 1 + nodes_per_layer + dum5 + dum2;
		    ++dum;
		}
		j = square_el_num - 3;
		*(connect + npel*dum)     = (dum3a+i)*square_el_num - i + j + dum5 + dum2;
		*(connect + npel*dum + 1) = (dum3a+i)*square_el_num - i + j + 1 + dum5 + dum2;
		*(connect + npel*dum + 2) = (dum3a+i+1)*square_el_num - i + j + dum5 + dum2;
		*(connect + npel*dum + 3) = (dum3a+i+1)*square_el_num - i + j - 1 + dum5 + dum2;
		*(connect + npel*dum + 4) = (dum3a+i)*square_el_num - i + j + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 5) = (dum3a+i)*square_el_num - i + j + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 6) = (dum3a+i+1)*square_el_num - i + j + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 7) = (dum3a+i+1)*square_el_num - i + j - 1 + nodes_per_layer + dum5 + dum2;
		++dum;
	    }
	    dum4 = square_el_num*dum3a;
	    *(connect + npel*dum)     = dum3a*(square_el_num-1) + square_el_num - 2 + dum5 + dum2;
	    *(connect + npel*dum + 1) = dum3a*(square_el_num-1) + square_el_num - 1 + dum5 + dum2;
	    *(connect + npel*dum + 2) = dum3a*(square_el_num-1) + square_el_num + dum5 + dum2;
	    *(connect + npel*dum + 3) = dum3a*square_el_num - 2 + square_el_num + dum5 + dum2;
	    *(connect + npel*dum + 4) = dum3a*(square_el_num-1) + square_el_num - 2 + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 5) = dum3a*(square_el_num-1) + square_el_num - 1 + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 6) = dum3a*(square_el_num-1) + square_el_num + nodes_per_layer + dum5 + dum2;
	    *(connect + npel*dum + 7) = dum3a*square_el_num - 2 + square_el_num + nodes_per_layer + dum5 + dum2;
	    ++dum;
	    for( i = 0; i < cylinder_el_num-1; ++i )
	    {
		*(connect + npel*dum)     = - i + (i+1)*square_el_num - 2 + dum4 + dum5 + dum2;
		*(connect + npel*dum + 1) = dum3a*(square_el_num-1) + square_el_num + i + dum5 + dum2;
		*(connect + npel*dum + 2) = dum3a*(square_el_num-1) + square_el_num + i + 1 + dum5 + dum2;
		*(connect + npel*dum + 3) = - i + (i+1)*square_el_num - 3 + dum4 + square_el_num + dum5 + dum2;
		*(connect + npel*dum + 4) = - i + (i+1)*square_el_num - 2 + dum4 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 5) = dum3a*(square_el_num-1) + square_el_num + i + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 6) = dum3a*(square_el_num-1) + square_el_num + i + 1 + nodes_per_layer + dum5 + dum2;
		*(connect + npel*dum + 7) = - i + (i+1)*square_el_num - 3 + dum4 + square_el_num + nodes_per_layer + dum5 + dum2;
		++dum;
	    }

	    dum2 += nodes_per_layer;
	}

/* These are the nodes for the 1st quarter */

	dum = 0;
	newnode_counter = 0;
	fdum3 = 0.0;
	fdum4 = square_length + square_el_length;
	for( k = 0; k < height_el_num+1; ++k )
	{
	    for( i = 0; i < square_el_num+1; ++i )
	    {
		for( j = 0; j < square_el_num+1; ++j )
		{
			fdum1 = square_el_length*((double)j);
			fdum2 = square_el_length*((double)i);
			*(coord + nsd*dum) = fdum1;
			*(coord + nsd*dum + 1) = fdum2;
			*(coord + nsd*dum + 2) = fdum3;
			++dum;
		
		}
		angle = ray_angle*((double)i);
		for( j = 0; j < cylinder_el_num; ++j )
		{
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
	    }
	    for( j = 0; j < cylinder_el_num; ++j )
	    {
		for( i = 0; i < square_el_num; ++i )
		{
		   angle = pi/2.0 - ray_angle*((double)i);
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
	    }

/* These are the nodes for the 2nd quarter */

	    for( i = 1; i < square_el_num+1; ++i )
	    {
		for( j = 0; j < square_el_num + 1; ++j )
		{
		   fdum1 = -square_el_length*((double)i);
		   fdum2 = square_el_length*((double)j);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
		angle = ray_angle*((double)i) + pi/2.0;
		for( j = 0; j < cylinder_el_num; ++j )
		{
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
	      }
	      for( j = 0; j < cylinder_el_num; ++j )
	      {
		for( i = 0; i < square_el_num; ++i )
		{
		   angle = pi/2.0 - ray_angle*((double)i) + pi/2.0;
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
	    }

/* These are the nodes for the 3rd quarter */

	    for( i = 1; i < square_el_num+1; ++i )
	    {
		for( j = 0; j < square_el_num + 1; ++j )
		{
		   fdum1 = -square_el_length*((double)j);
		   fdum2 = -square_el_length*((double)i);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
		angle = ray_angle*((double)i) + pi;
		for( j = 0; j < cylinder_el_num; ++j )
		{
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
	      }
	      for( j = 0; j < cylinder_el_num; ++j )
	      {
		for( i = 0; i < square_el_num; ++i )
		{
		   angle = pi/2.0 - ray_angle*((double)i) + pi;
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
	    }

/* These are the nodes for the 4th quarter */

	    for( i = 1; i < square_el_num+1; ++i )
	    {
		for( j = 1; j < square_el_num + 1; ++j )
		{
		   fdum1 = square_el_length*((double)i);
		   fdum2 = -square_el_length*((double)j);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
		angle = ray_angle*((double)i) + 3.0*pi/2.0;
		for( j = 0; j < cylinder_el_num; ++j )
		{
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		}
	    }
	    j = 0;
	    for( i = 1; i < square_el_num; ++i )
	    {
		angle = pi/2.0 - ray_angle*((double)i) + 3.0*pi/2.0;
		fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		*(coord + nsd*dum) = fdum1;
		*(coord + nsd*dum + 1) = fdum2;
		*(coord + nsd*dum + 2) = fdum3;
		++dum;
	    }
	    for( j = 1; j < cylinder_el_num; ++j )
	    {
		for( i = 1; i < square_el_num-1; ++i )
		{
		   angle = pi/2.0 - ray_angle*((double)i) + 3.0*pi/2.0;
		   fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		   fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		   *(coord + nsd*dum) = fdum1;
		   *(coord + nsd*dum + 1) = fdum2;
		   *(coord + nsd*dum + 2) = fdum3;
		   ++dum;
		   if( k == half_height_el_num )
		   {
/* These are the new inserted nodes.  They are inserted at the middle of the mesh
   to create a nick.  They are added immediately on top of the previous node
   at the same location.
*/

			angle = pi/2.0 - ray_angle*((double)i) + 3.0*pi/2.0;
			fdum1 =
			    (ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
			fdum2 =
			    (ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
			/*fprintf( o1, "%5d  ",dum);
			fprintf( o1, "jjjj %14.6e  %14.6e  %14.6e",fdum1,fdum2,fdum3);
			fprintf( o1, "\n");*/
			*(coord + nsd*dum) = fdum1;
			*(coord + nsd*dum + 1) = fdum2;
			*(coord + nsd*dum + 2) = fdum3;

/* I record the new nodes which are necessary to modify the connectivity. */

			newnode[newnode_counter].node = dum;
			++dum;
			++newnode_counter;
		   }
		}
		i = square_el_num-1;
		angle = pi/2.0 - ray_angle*((double)i) + 3.0*pi/2.0;
		fdum1 =
			(ratio*fdum4+cylinder_el_length*((double)j))*cos(angle);
		fdum2 =
			(ratio*fdum4+cylinder_el_length*((double)j))*sin(angle);
		*(coord + nsd*dum) = fdum1;
		*(coord + nsd*dum + 1) = fdum2;
		*(coord + nsd*dum + 2) = fdum3;
		++dum;
	    }

	    fdum3 += height_el_length;
	    
	}

	printf( "%6dzzzzzz\n", dum);

/* Re-order the newnodes into accending order.  This is needed based on the way
   the connectivity is adjusted upward for the inserted nodes.
*/

	for( i = 0; i < newnode_counter; ++i ) newnode[i].index = i;
	qsort( newnode, newnode_counter, sizeof(NEWNODE), compare_newnode);
	for( i = 0; i < newnode_counter; ++i )
	{
		printf( "%6d %9d\n", newnode[i].index, newnode[i].node);
	}


/* Shift the connectivity to reflect the inserted nodes. */

	for( i = 0; i < numel; ++i )
	{
	    for( j = 0; j < npel; ++j )
	    {
		for( k = 0; k < newnode_counter; ++k )
		{
		    if(*(connect + npel*i + j) >= newnode[k].node) *(connect + npel*i + j) += 1;
		    if( i > half_height_el_num*nodes_per_layer) 
		    {
/* This code is for re-doing the connectivity of the elements that have the inserted nodes
   themselves.  The reason it is done for elements that are

      i > half_height_el_num*nodes_per_layer

   (i.e. elements which are above the middle height) is because the element layer below
   should continue to use the non-inserted nodes.  This is a significant consideration
   in the future since I have to decide which element keeps the old node and which gets
   the new node inserted on top.  I have thought about many ways to do this, like comparing
   adjacent elements which share inserted nodes vs. those that do not, and assigning based
   on that criteria.

   Other factors may be that I will insert nodes on an element by element basis, such as
   when the stress at a node goes beyond a certain point.  This may make choosing how
   to associate the new nodes easier.
*/

   
			if(*(connect + npel*i + j) == newnode[k].node - 1)
			    *(connect + npel*i + j) = newnode[k].node;
		    }
		}
	    }
	}

	fprintf( o1, "el no., connectivity, matl no. \n");
	for( i = 0; i < numel; ++i )
	{
	    fprintf( o1, "%6d %9d %9d %9d %9d %9d %9d %9d %9d %9d", i,
		    *(connect + npel*i),
		    *(connect + npel*i + 1),
		    *(connect + npel*i + 2),
		    *(connect + npel*i + 3),
		    *(connect + npel*i + 4),
		    *(connect + npel*i + 5),
		    *(connect + npel*i + 6),
		    *(connect + npel*i + 7),
		    0);
	    fprintf( o1, "\n");
	    ++dum;
	}

	fprintf( o1, "node no., coordinates \n");
	for( k = 0; k < numnp; ++k )
	{
		fprintf( o1, "%5d  ",k);
		fprintf( o1, "%14.6e  %14.6e  %14.6e",
			*(coord + nsd*k), *(coord + nsd*k + 1), *(coord + nsd*k + 2));
		fprintf( o1, "\n");
	}

        dum= 0;
        dum2 = nodes_per_layer*height_el_num + (cylinder_el_num)*(square_el_num-2);
        fprintf( o1, "prescribed displacement x: node  disp value\n");
        for( k = dum2; k < numnp; ++k )
	{
                fprintf( o1, "%4d %12.6e\n",k,0.0);
        }
        fprintf( o1, "%4d\n ",-10);
        fprintf( o1, "prescribed displacement y: node  disp value\n");
        for( k = dum2; k < numnp; ++k )
	{
                fprintf( o1, "%4d %12.6e\n",k,0.0);
        }
        fprintf( o1, "%4d\n ",-10);
        fprintf( o1, "prescribed displacement z: node  disp value\n");
        for( k = dum2; k < numnp; ++k )
	{
                fprintf( o1, "%4d %12.6e\n",k,0.0);
        }
        fprintf( o1, "%4d\n ",-10);

	dum2 = nodes_per_layer;
	dum3 = (square_el_num+cylinder_el_num+1);
	dum3a = (square_el_num+cylinder_el_num);
	fdum4 = 2.0*square_length*square_length;
	fdum4 = sqrt(fdum4);

        fprintf( o1, "node with point load and load vector in x,y\n");
	for( i = 0; i < nodes_per_layer; ++i )
	{
		fprintf( o1, "%4d %14.6e %14.6e %14.6e\n", i, 0.0, 0.0, -1.0e9);
	}


        fprintf( o1, "%4d\n ",-10);
        fprintf( o1, "node no. with stress and stress vector in xx,yy,xy,zx,yz\n");
        fprintf( o1, "%4d ",-10);

        return 1;
}

