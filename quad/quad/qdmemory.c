/*
    This utility function allocates the memory for
    the quad finite element program.

        Updated 6/4/00

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
#include "qdconst.h"
#include "qdstruct.h"

int qdMemory( double **mem_double, int sofmf, int **mem_int, int sofmi,
	MATL **matl, int nmat, XYI **mem_XYI, int sofmXYI, SDIM **mem_SDIM,
	int sofmSDIM, STRAIN **strain, STRESS **stress, int sofmSTRESS )
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
/* For the STRAIN doubles */

	*strain=(STRAIN *)calloc(sofmSTRESS,sizeof(STRAIN));
	if(!strain )
	{
		printf( "failed to allocate memory for strain doubles\n ");
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

/* For the SDIM doubles */

	*mem_SDIM=(SDIM *)calloc(sofmSDIM,sizeof(SDIM));
	if(!mem_SDIM )
	{
	        printf( "failed to allocate memory for SDIM doubles\n ");
	        exit(1);
	}

        return 1;
}


int qdReGetMemory( double **mem_double, int sofmf, int **mem_int, int sofmi,
	MATL **matl, int nmat, XYI **mem_XYI, int sofmXYI, SDIM **mem_SDIM,
	int sofmSDIM, STRAIN **strain, STRESS **stress, int sofmSTRESS )
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

/* For the STRAIN doubles */

	*strain=(STRAIN *)realloc(*strain, sofmSTRESS*sizeof(STRAIN));
	if(!strain )
	{
		printf( "failed to allocate memory for strain doubles\n ");
		exit(1);
	}
	memset(*strain,0,sofmSTRESS*sizeof(STRAIN));

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

/* For the SDIM doubles */

	*mem_SDIM=(SDIM *)realloc(*mem_SDIM, sofmSDIM*sizeof(SDIM));
	if(!mem_SDIM )
	{
	        printf( "failed to allocate memory for SDIM doubles\n ");
	        exit(1);
	}
	memset(*mem_SDIM,0,sofmSDIM*sizeof(SDIM));

	return 1;
}

