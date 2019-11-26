/*
    This library function reads in additional stress and strain
    data for the graphics program for wedge elements.

		Updated 3/19/01

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "../wedge/weconst.h"
#if WEDGE1
#include "../wedge/westruct.h"
#endif
#if WEDGE2
#include "../wedge2/we2struct.h"
#endif

extern int dof, nmat, numel, numnp;

int wereader_gr( FILE *o1, SDIM *strain_node, SDIM *stress_node)
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

