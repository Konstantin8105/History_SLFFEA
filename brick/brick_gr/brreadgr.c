/*
    This library function reads in additional stress and strain
    data for the graphics program for brick elements.

		Updated 12/14/00

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "../brick/brconst.h"
#if BRICK1
#include "../brick/brstruct.h"
#endif
#if BRICK2
#include "../brick2/br2struct.h"
#endif

extern int dof, nmat, numel, numnp;

int brreader_gr( FILE *o1, SDIM *strain_node, SDIM *stress_node)
{
        int i,j,dum;
	char buf[ BUFSIZ ];
	char text;

        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );

        printf("principal stress for node: ");
        fscanf( o1,"%d",&dum);
        printf( "(%4d)",dum);
        while( dum > -1 )
        {
           fscanf( o1,"%lf ",&stress_node[dum].I);
           fscanf( o1,"%lf ",&stress_node[dum].II);
           fscanf( o1,"%lf ",&stress_node[dum].III);
           printf(" %12.5e",stress_node[dum].I);
           printf(" %12.5e",stress_node[dum].II);
           printf(" %12.5e",stress_node[dum].III);
           fscanf( o1,"\n");
           printf( "\n");
           printf("principal stress for node: ");
           fscanf( o1,"%d",&dum);
           printf( "(%4d)",dum);
        }
        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );
        printf( "\n\n");

        printf("strain for node: ");
        fscanf( o1,"%d",&dum);
        printf( "(%4d)",dum);
        while( dum > -1 )
        {
           fscanf( o1,"%lf ",&strain_node[dum].xx);
           fscanf( o1,"%lf ",&strain_node[dum].yy);
           fscanf( o1,"%lf ",&strain_node[dum].zz);
           fscanf( o1,"%lf ",&strain_node[dum].xy);
           fscanf( o1,"%lf ",&strain_node[dum].zx);
           fscanf( o1,"%lf ",&strain_node[dum].yz);
           printf(" %12.5e",strain_node[dum].xx);
           printf(" %12.5e",strain_node[dum].yy);
           printf(" %12.5e",strain_node[dum].zz);
           printf(" %12.5e",strain_node[dum].xy);
           printf(" %12.5e",strain_node[dum].zx);
           printf(" %12.5e",strain_node[dum].yz);
           fscanf( o1,"\n");
           printf( "\n");
           printf("strain for node: ");
           fscanf( o1,"%d",&dum);
           printf( "(%4d)",dum);
        }
        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );
        printf( "\n\n");

        printf("principal strain for node: ");
        fscanf( o1,"%d",&dum);
        printf( "(%4d)",dum);
        while( dum > -1 )
        {
           fscanf( o1,"%lf ",&strain_node[dum].I);
           fscanf( o1,"%lf ",&strain_node[dum].II);
           fscanf( o1,"%lf ",&strain_node[dum].III);
           printf(" %12.5e",strain_node[dum].I);
           printf(" %12.5e",strain_node[dum].II);
           printf(" %12.5e",strain_node[dum].III);
           fscanf( o1,"\n");
           printf( "\n");
           printf("principal strain for node: ");
           fscanf( o1,"%d",&dum);
           printf( "(%4d)",dum);
        }
        printf( "\n\n");

        return 1;
}

