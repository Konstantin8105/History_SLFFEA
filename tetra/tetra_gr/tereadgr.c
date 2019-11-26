/*
    This library function reads in additional stress and strain
    data for the graphics program for tetrahedral elements.

		Updated 9/5/01

    SLFFEA source file
    Version:  1.3
    Copyright (C) 1999, 2000, 2001, 2002  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "../tetra/teconst.h"
#if TETRA1
#include "../tetra/testruct.h"
#endif
#if TETRA2
#include "../tetra2/te2struct.h"
#endif

extern int dof, nmat, numel, numnp;

int tereader_gr( FILE *o1, STRAIN *strain_node, STRESS *stress_node)
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

