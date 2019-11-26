/*
    This library function writes the connectivity data for surface
    elements for a finite element program which does analysis on a brick 

		Updated 11/11/00

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "brconst.h"

extern int dof, nmat, nmode, numel, numnp, gauss_stress_print_flag;
extern int numel_surf;

int brConnectSurfwriter ( int *connect_surf, int *el_matl_surf, char *name)
{
        int i,j,dum,check, name_length;
	char *ccheck;
	char connect_dat[30];
	FILE *o5;

/* Every output file is named after the input file with
   a ".obr" extension */

	name_length = strlen(name);
	if( name_length > 25) name_length = 25;

/* Open surface connectivity data file */

	memset(connect_dat,0,30*sizeof(char));

	ccheck = strncpy(connect_dat, name, name_length);
	if(!ccheck) printf( " Problems with strncpy \n");

	ccheck = strncpy(connect_dat+name_length, ".con", 4);
	if(!ccheck) printf( " Problems with strncpy \n");

	o5 = fopen( connect_dat,"w" );

/* print out the connectivity for the surface elements */

	fprintf( o5, "   numel_surf  Surface el. conectivity(This is for the brick file: %s) \n ", name);
	fprintf( o5, "    %4d \n ", numel_surf);
	fprintf( o5, "el no., connectivity, matl no. \n");

	for( i = 0; i < numel_surf; ++i )
	{
		fprintf( o5, "%6d ",i);
		for( j = 0; j < npel; ++j )
		{
			fprintf( o5, "%6d ",*(connect_surf+npel*i+j));
		}
		fprintf( o5, "   %3d\n",*(el_matl_surf+i));
	}

        return 1;
}

