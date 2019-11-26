/*
    This library function reads in additional stress and strain
    data for the graphics program for elements with additional
    stresses I, II.

		Updated 8/12/06

    SLFFEA source file
    Version:  1.4
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#if QUAD1
#include "../quad/quad/qdconst.h"
#include "../quad/quad/qdstruct.h"
#endif
#if TRI1
#include "../tri/tri/trconst.h"
#include "../tri/tri/trstruct.h"
#endif

extern int dof, nmat, numel, numnp;

int reader_gr2( FILE *o1, SDIM *strain_node, SDIM *stress_node)
{
	int i,j,dum;
	char buf[ BUFSIZ ];
	char text;

	fscanf( o1,"\n");
	fgets( buf, BUFSIZ, o1 );

	printf("principal stress for node: ");
	fscanf( o1,"%d",&dum);
	printf( "(%6d)",dum);
	while( dum > -1 )
	{
	   fscanf( o1,"%lf ",&stress_node[dum].I);
	   fscanf( o1,"%lf ",&stress_node[dum].II);
	   printf(" %12.5e",stress_node[dum].I);
	   printf(" %12.5e",stress_node[dum].II);
	   fscanf( o1,"\n");
	   printf( "\n");
	   printf("principal stress for node: ");
	   fscanf( o1,"%d",&dum);
	   printf( "(%6d)",dum);
	}
	fscanf( o1,"\n");
	fgets( buf, BUFSIZ, o1 );
	printf( "\n\n");

	printf("strain for node: ");
	fscanf( o1,"%d",&dum);
	printf( "(%6d)",dum);
	while( dum > -1 )
	{
	   fscanf( o1,"%lf ",&strain_node[dum].xx);
	   fscanf( o1,"%lf ",&strain_node[dum].yy);
	   fscanf( o1,"%lf ",&strain_node[dum].xy);
	   printf(" %12.5e",strain_node[dum].xx);
	   printf(" %12.5e",strain_node[dum].yy);
	   printf(" %12.5e",strain_node[dum].xy);
	   fscanf( o1,"\n");
	   printf( "\n");
	   printf("strain for node: ");
	   fscanf( o1,"%d",&dum);
	   printf( "(%6d)",dum);
	}
	fscanf( o1,"\n");
	fgets( buf, BUFSIZ, o1 );
	printf( "\n\n");

	printf("principal strain for node: ");
	fscanf( o1,"%d",&dum);
	printf( "(%6d)",dum);
	while( dum > -1 )
	{
	   fscanf( o1,"%lf ",&strain_node[dum].I);
	   fscanf( o1,"%lf ",&strain_node[dum].II);
	   printf(" %12.5e",strain_node[dum].I);
	   printf(" %12.5e",strain_node[dum].II);
	   fscanf( o1,"\n");
	   printf( "\n");
	   printf("principal strain for node: ");
	   fscanf( o1,"%d",&dum);
	   printf( "(%6d)",dum);
	}
	printf( "\n\n");

	return 1;
}
