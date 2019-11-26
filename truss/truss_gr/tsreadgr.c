/*
    This library function reads in additional strain
    data for the graphics program for truss elements.

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
#include "../truss/tsconst.h"
#include "../truss/tsstruct.h"

extern int dof, nmat, numel, numnp;

int tsreader_gr( FILE *o1, STRAIN *strain)
{
        int i,j,dum;
	char buf[ BUFSIZ ];
	char text;

        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );

        printf("strain for ele: ");
        fscanf( o1,"%d",&dum);
        printf( "(%4d)",dum);
        while( dum > -1 )
        {
           fscanf( o1,"%lf ",&strain[dum].xx);
           printf(" %12.5e",strain[dum].xx);
           fscanf( o1,"\n");
           printf( "\n");
           printf("strain for ele: ");
           fscanf( o1,"%d",&dum);
           printf( "(%4d)",dum);
        }
        printf( "\n\n");

        return 1;
}

