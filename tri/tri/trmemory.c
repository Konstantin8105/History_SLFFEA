/*
    This utility function allocates the memory for
    the triangle finite element program.

        Updated 8/13/01

    SLFFEA source file
    Version:  1.3
    Copyright (C) 1999, 2000, 2001, 2002  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/


#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "trconst.h"
#include "trstruct.h"

int trMemory( double **mem_double, int sofmf, int **mem_int, int sofmi,
	MATL **matl, int nmat, XYI **mem_XYI, int sofmXYI, STRAIN **strain,
	STRAIN **strain_node, STRESS **stress, STRESS **stress_node,
	int sofmSTRESS, int sofmSTRESS_node )
{
/* For the doubles */
	*mem_double=(double *)calloc(sofmf,sizeof(double));

	if(!mem_double )
	{
	    printf( "failed to allocate memory for doubles\n ");
	    exit(1);
	}

/* For the materials */
	*matl=(MATL *)calloc(nmat,sizeof(MATL));
	if(!matl )
	{
		printf( "failed to allocate memory for matl doubles\n ");
		exit(1);
	}

/* For the STRESS doubles */

	*stress=(STRESS *)calloc(sofmSTRESS,sizeof(STRESS));
	if(!stress )
	{
		printf( "failed to allocate memory for stress doubles\n ");
		exit(1);
	}

/* For the STRESS node doubles */

	*stress_node=(STRESS *)calloc(sofmSTRESS_node,sizeof(STRESS));
	if(!stress_node )
	{
		printf( "failed to allocate memory for stress_node doubles\n ");
		exit(1);
	}

/* For the STRAIN doubles */

	*strain=(STRAIN *)calloc(sofmSTRESS,sizeof(STRAIN));
	if(!strain )
	{
		printf( "failed to allocate memory for strain doubles\n ");
		exit(1);
	}

/* For the STRAIN node doubles */

	*strain_node=(STRAIN *)calloc(sofmSTRESS_node,sizeof(STRAIN));
	if(!strain_node )
	{
		printf( "failed to allocate memory for strain_node doubles\n ");
		exit(1);
	}

/* For the integers */
	*mem_int=(int *)calloc(sofmi,sizeof(int));
	if(!mem_int )
	{
		printf( "failed to allocate memory for integers\n ");
		exit(1);
	}

/* For the XYI integers */

	*mem_XYI=(XYI *)calloc(sofmXYI,sizeof(XYI));
	if(!mem_XYI )
	{
		printf( "failed to allocate memory for XYI integers\n ");
		exit(1);
	}

        return 1;
}


int trReGetMemory( double **mem_double, int sofmf, int **mem_int, int sofmi,
	MATL **matl, int nmat, XYI **mem_XYI, int sofmXYI, STRAIN **strain,
	STRAIN **strain_node, STRESS **stress, STRESS **stress_node,
	int sofmSTRESS, int sofmSTRESS_node )
{
/* For the doubles */
	*mem_double=(double *)realloc(*mem_double, sofmf*sizeof(double));

	if(!mem_double )
	{
	    printf( "failed to allocate memory for doubles\n ");
	    exit(1);
	}
	memset(*mem_double,0,sofmf*sizeof(double));

/* For the materials */
	*matl=(MATL *)realloc(*matl, nmat*sizeof(MATL));
	if(!matl )
	{
		printf( "failed to allocate memory for matl doubles\n ");
		exit(1);
	}
	memset(*matl,0,nmat*sizeof(MATL));

/* For the STRESS doubles */

	*stress=(STRESS *)realloc(*stress, sofmSTRESS*sizeof(STRESS));
	if(!stress )
	{
		printf( "failed to allocate memory for stress doubles\n ");
		exit(1);
	}
	memset(*stress,0,sofmSTRESS*sizeof(STRESS));

/* For the STRESS node doubles */

	*stress_node=(STRESS *)realloc(*stress_node, sofmSTRESS_node*sizeof(STRESS));
	if(!stress_node )
	{
		printf( "failed to allocate memory for stress_node doubles\n ");
		exit(1);
	}
	memset(*stress_node,0,sofmSTRESS_node*sizeof(STRESS));

/* For the STRAIN doubles */

	*strain=(STRAIN *)realloc(*strain, sofmSTRESS*sizeof(STRAIN));
	if(!strain )
	{
		printf( "failed to allocate memory for strain doubles\n ");
		exit(1);
	}
	memset(*strain,0,sofmSTRESS*sizeof(STRAIN));

/* For the STRAIN node doubles */

	*strain_node=(STRAIN *)realloc(*strain_node, sofmSTRESS_node*sizeof(STRAIN));
	if(!strain_node )
	{
		printf( "failed to allocate memory for strain_node doubles\n ");
		exit(1);
	}
        memset(*strain_node,0,sofmSTRESS_node*sizeof(STRAIN));

/* For the integers */
	*mem_int=(int *)realloc(*mem_int, sofmi*sizeof(int));
	if(!mem_int )
	{
		printf( "failed to allocate memory for integers\n ");
		exit(1);
	}
	memset(*mem_int,0,sofmi*sizeof(int));

/* For the XYI integers */

	*mem_XYI=(XYI *)realloc(*mem_XYI, sofmXYI*sizeof(XYI));
	if(!mem_XYI )
	{
		printf( "failed to allocate memory for XYI integers\n ");
		exit(1);
	}
	memset(*mem_XYI,0,sofmXYI*sizeof(XYI));

	return 1;
}

