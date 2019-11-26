/*
    This library function reads in data for a finite element
    program which does analysis on a plate element

		Updated 10/17/06

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
#include "plconst.h"
#include "plstruct.h"

extern int dof, nmat, nmode, numel, numnp;
extern int stress_read_flag, element_stress_read_flag;

int plreader( BOUND bc, int *connect, double *coord, int *el_matl, double *force,
	MATL *matl, MOMENT *moment, MDIM *moment_node, char *name, FILE *o1,
	STRESS *stress, SDIM *stress_node, double *U)
{
        int i,j,dum,dum2, name_length;
	char *ccheck;
	char buf[ BUFSIZ ];
	char text, stress_dat[30];
	FILE *o4;

	if(element_stress_read_flag)
	{
/* Open stress data output file */

		name_length = strlen(name);
		if( name_length > 25) name_length = 25;

		memset(stress_dat,0,30*sizeof(char));

		ccheck = strncpy(stress_dat, name, name_length);
		if(!ccheck) printf( " Problems with strncpy \n");

		ccheck = strncpy(stress_dat+name_length, ".str", 4);
		if(!ccheck) printf( " Problems with strncpy \n");

		o4 = fopen( stress_dat,"r" );
		if(o4 == NULL ) {
			printf("Can't find file %30s\n",stress_dat);
			element_stress_read_flag = 0;
		}
	}

        printf( "number of elements:%d nodes:%d materials:%d modes:%d dof:%d\n",
                numel,numnp,nmat,nmode,dof);
        fgets( buf, BUFSIZ, o1 );
        printf( "\n");

        for( i = 0; i < nmat; ++i )
        {
           fscanf( o1, "%d ",&dum);
           printf( "material (%3d) Emod, nu, density, thickness, shear fac.: ",dum);
           fscanf( o1, " %lf %lf %lf %lf %lf\n", &matl[dum].E, &matl[dum].nu,
		&matl[dum].rho, &matl[dum].thick, &matl[dum].shear);
           printf( " %9.4f %9.4f %9.4f %9.4f %9.4f\n", matl[dum].E, matl[dum].nu,
		matl[dum].rho, matl[dum].thick, matl[dum].shear);
        }
        fgets( buf, BUFSIZ, o1 );
        printf( "\n");

        for( i = 0; i < numel; ++i )
        {
           fscanf( o1,"%d ",&dum);
           printf( "connectivity for element (%4d) ",dum);
           for( j = 0; j < npel; ++j )
           {
                fscanf( o1, "%d",(connect+npel*dum+j));
                printf( "%4d ",*(connect+npel*dum+j));
           }
           fscanf( o1,"%d\n",(el_matl+dum));
           printf( " with matl %3d\n",*(el_matl+dum));
        }
        fgets( buf, BUFSIZ, o1 );
        printf( "\n");

        for( i = 0; i < numnp; ++i )
        {
           fscanf( o1,"%d ",&dum);
           printf( "coordinate (%d) ",dum);
           printf( "coordinates ");
           for( j = 0; j < nsd; ++j )
           {
                fscanf( o1, "%lf ",(coord+nsd*dum+j));
                printf( "%9.4f ",*(coord+nsd*dum+j));
           }
           fscanf( o1,"\n");
           printf( "\n");
        }
        fgets( buf, BUFSIZ, o1 );
        printf( "\n");

        dum= 0;
        fscanf( o1,"%d",&bc.fix[dum].z);
        printf( "node (%4d) has a z prescribed displacement of: ",bc.fix[dum].z);
        while( bc.fix[dum].z > -1 )
        {
                fscanf( o1,"%lf\n%d",(U+ndof*bc.fix[dum].z),
                        &bc.fix[dum+1].z);
                printf( "%14.6e\n",*(U+ndof*bc.fix[dum].z));
                printf( "node (%4d) has a z prescribed displacement of: ",
                        bc.fix[dum+1].z);
                ++dum;
        }
        bc.num_fix[0].z=dum;
        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );
        printf( "\n\n");

        dum= 0;
        fscanf( o1,"%d",&bc.fix[dum].phix);
        printf( "node (%4d) has a prescribed angle phi x of: ",bc.fix[dum].phix);
        while( bc.fix[dum].phix > -1 )
        {
                fscanf( o1,"%lf\n%d",(U+ndof*bc.fix[dum].phix+1),
                        &bc.fix[dum+1].phix);
                printf( "%14.6e\n",*(U+ndof*bc.fix[dum].phix+1));
                printf( "node (%4d) has a prescribed angle phi x of: ",
                        bc.fix[dum+1].phix);
                ++dum;
        }
        bc.num_fix[0].phix=dum;
        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );
        printf( "\n\n");

        dum= 0;
        fscanf( o1,"%d",&bc.fix[dum].phiy);
        printf( "node (%4d) has a prescribed angle phi y of: ",bc.fix[dum].phiy);
        while( bc.fix[dum].phiy > -1 )
        {
                fscanf( o1,"%lf\n%d",(U+ndof*bc.fix[dum].phiy+2),
                        &bc.fix[dum+1].phiy);
                printf( "%14.6e\n",*(U+ndof*bc.fix[dum].phiy+2));
                printf( "node (%4d) has a prescribed angle phi y of: ",
                        bc.fix[dum+1].phiy);
                ++dum;
        }
        bc.num_fix[0].phiy=dum;
        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );
        printf( "\n\n");

        dum= 0;
        printf("force vector for node: ");
        fscanf( o1,"%d",&bc.force[dum]);
        printf( "(%4d)",bc.force[dum]);
        while( bc.force[dum] > -1 )
        {
           for( j = 0; j < ndof; ++j )
           {
                fscanf( o1,"%lf ",(force+ndof*bc.force[dum]+j));
                printf("%16.4f ",*(force+ndof*bc.force[dum]+j));
           }
           fscanf( o1,"\n");
           printf( "\n");
           printf("force vector for node: ");
           ++dum;
           fscanf( o1,"%d",&bc.force[dum]);
           printf( "(%4d)",bc.force[dum]);
        }
        bc.num_force[0]=dum;
        fscanf( o1,"\n");
        fgets( buf, BUFSIZ, o1 );
        printf( "\n\n");

	if(stress_read_flag)
	{
	   printf("stress for node: ");
	   fscanf( o1,"%d",&dum);
	   printf( "(%4d)",dum);
	   while( dum > -1 )
	   {
		fscanf( o1,"%lf ",&moment_node[dum].xx);
		fscanf( o1,"%lf ",&moment_node[dum].yy);
		fscanf( o1,"%lf ",&moment_node[dum].xy);
		fscanf( o1,"%lf ",&stress_node[dum].zx);
		fscanf( o1,"%lf ",&stress_node[dum].yz);
		printf(" %12.5e",moment_node[dum].xx);
		printf(" %12.5e",moment_node[dum].yy);
		printf(" %12.5e",moment_node[dum].xy);
		printf(" %12.5e",stress_node[dum].zx);
		printf(" %12.5e",stress_node[dum].yz);
		fscanf( o1,"\n");
		printf( "\n");
		printf("stress for node: ");
		fscanf( o1,"%d",&dum);
		printf( "(%4d)",dum);
	   }
	}
	printf( "\n\n");

	if(element_stress_read_flag)
        {
	   fgets( buf, BUFSIZ, o4 );
	   printf( "\n\n");
	   printf("stress for ele: ");
	   fscanf( o4,"%d",&dum);
	   printf( "(%4d)",dum);
           while( dum > -1 )
           {
		fscanf( o4,"%d",&dum2);
		printf( " node (%1d)",dum2);
		fscanf( o4,"%lf ",&moment[dum].pt[dum2].xx);
		fscanf( o4,"%lf ",&moment[dum].pt[dum2].yy);
		fscanf( o4,"%lf ",&moment[dum].pt[dum2].xy);
		fscanf( o4,"%lf ",&stress[dum].pt[dum2].zx);
		fscanf( o4,"%lf ",&stress[dum].pt[dum2].yz);
		printf(" %12.5e",moment[dum].pt[dum2].xx);
		printf(" %12.5e",moment[dum].pt[dum2].yy);
		printf(" %12.5e",moment[dum].pt[dum2].xy);
		printf(" %12.5e",stress[dum].pt[dum2].zx);
		printf(" %12.5e",stress[dum].pt[dum2].yz);
		fscanf( o4,"\n");
		printf( "\n");
		printf("stress for ele: ");
		fscanf( o4,"%d",&dum);
		printf( "(%4d)",dum);
           }
        }
        printf( "\n\n");

        return 1;
}

