/*
    This library function reads in additional strain and
    curvature data for the graphics program for beam elements.

		Updated 12/14/99

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "../beam/bmconst.h"
#include "../beam/bmstruct.h"

extern int dof, nmat, numel, numnp;

int bmreader_gr( FILE *o1, CURVATURE *curve, STRAIN *strain)
{
        int i,j,dum,dum2;
	char buf[ BUFSIZ ];
	char text;

        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );

        printf("strain for ele: ");
        fscanf( o1,"%d",&dum);
        printf( "(%4d)",dum);
        while( dum > -1 )
        {
           fscanf( o1,"%d",&dum2);
           printf( " node (%1d)",dum2);
           fscanf( o1,"%lf ",&strain[dum].pt[dum2].xx);
           fscanf( o1,"%lf ",&curve[dum].pt[dum2].xx);
           fscanf( o1,"%lf ",&curve[dum].pt[dum2].yy);
           fscanf( o1,"%lf ",&curve[dum].pt[dum2].zz);
           printf(" %12.5e",strain[dum].pt[dum2].xx);
           printf(" %12.5e",curve[dum].pt[dum2].xx);
           printf(" %12.5e",curve[dum].pt[dum2].yy);
           printf(" %12.5e",curve[dum].pt[dum2].zz);
           fscanf( o1,"\n");
           printf( "\n");
           printf("strain for ele: ");
           fscanf( o1,"%d",&dum);
           printf( "(%4d)",dum);
        }
        printf( "\n\n");

        return 1;
}

