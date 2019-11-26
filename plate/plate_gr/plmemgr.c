/*
    This utility function allocates the memory for
    the plate finite element graphics program.

  
   			Last Update 3/16/00

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
#include "../plate/plconst.h"
#include "../plate/plstruct.h"
#include "plstrcgr.h"

extern int input_flag, post_flag;

int plMemory_gr( ICURVATURE **curve_color, IMOMENT **moment_color,
	ISTRAIN **strain_color, ISTRESS **stress_color, NORM **mem_NORM,
	int sofmISTRESS,  int sofmNORM )
{
/* Allocate the memory for when first running the program

        Updated 12/23/99

*/

/* For the IMOMENT integers */

        *moment_color=(IMOMENT *)calloc(sofmISTRESS,sizeof(IMOMENT));
        if(!moment_color )
        {
                printf( "failed to allocate memory for moment integers\n ");
                exit(1);
        }

/* For the ISTRESS integers */

        *stress_color=(ISTRESS *)calloc(sofmISTRESS,sizeof(ISTRESS));
        if(!stress_color )
        {
                printf( "failed to allocate memory for stress integers\n ");
                exit(1);
        }

/* For the ICURVATURE integers */

        *curve_color=(ICURVATURE *)calloc(sofmISTRESS,sizeof(ICURVATURE));
        if(!curve_color )
        {
                printf( "failed to allocate memory for curve integers\n ");
                exit(1);
        }

/* For the ISTRAIN integers */

        *strain_color=(ISTRAIN *)calloc(sofmISTRESS,sizeof(ISTRAIN));
        if(!strain_color )
        {
                printf( "failed to allocate memory for strain integers\n ");
                exit(1);
        }

/* For the NORM doubles */

	*mem_NORM=(NORM *)calloc(sofmNORM, sizeof(NORM));
	if(!mem_NORM )
	{
		printf( "failed to allocate memory for NORM doubles\n ");
		exit(1);
	}

  	return 1; 
}


int plMemory2_gr( ZPhiF **mem_ZPhiF, int sofmZPhiF)
{

/* For the ZPhiF doubles */

        *mem_ZPhiF=(ZPhiF *)calloc(sofmZPhiF,sizeof(ZPhiF));
        if(!mem_ZPhiF )
        {
                printf( "failed to allocate memory for ZPhiF doubles\n ");
                exit(1);
        }

  	return 1; 
}

int plReGetMemory_gr( ICURVATURE **curve_color, IMOMENT **moment_color,
	ISTRAIN **strain_color, ISTRESS **stress_color, NORM **mem_NORM,
	int sofmISTRESS, int sofmNORM )
{
/* Allocate the memory for when first running the program

        Updated 3/14/00

*/

/* For the IMOMENT integers */

        *moment_color=(IMOMENT *)realloc(*moment_color, sofmISTRESS*sizeof(IMOMENT));
        if(!moment_color )
        {
                printf( "failed to allocate memory for moment integers\n ");
                exit(1);
        }
	memset(*moment_color,0,sofmISTRESS*sizeof(IMOMENT));
	

/* For the ISTRESS integers */

        *stress_color=(ISTRESS *)realloc(*stress_color, sofmISTRESS*sizeof(ISTRESS));
        if(!stress_color )
        {
                printf( "failed to allocate memory for stress integers\n ");
                exit(1);
        }
	memset(*stress_color,0,sofmISTRESS*sizeof(ISTRESS));

/* For the ICURVATURE integers */

        *curve_color=(ICURVATURE *)realloc(*curve_color, sofmISTRESS*sizeof(ICURVATURE));
        if(!curve_color )
        {
                printf( "failed to allocate memory for curve integers\n ");
                exit(1);
        }
	memset(*curve_color,0,sofmISTRESS*sizeof(ICURVATURE));

/* For the ISTRAIN integers */

        *strain_color=(ISTRAIN *)realloc(*strain_color, sofmISTRESS*sizeof(ISTRAIN));
        if(!strain_color )
        {
                printf( "failed to allocate memory for strain integers\n ");
                exit(1);
        }
	memset(*strain_color,0,sofmISTRESS*sizeof(ISTRAIN));

/* For the NORM doubles */

	*mem_NORM=(NORM *)realloc(*mem_NORM, sofmNORM*sizeof(NORM));
	if(!mem_NORM )
	{
                printf( "failed to allocate memory for NORM doubles\n ");
                exit(1);
	}
	memset(*mem_NORM,0,sofmNORM*sizeof(NORM));

  	return 1; 
}


int plReGetMemory2_gr( ZPhiF **mem_ZPhiF, int sofmZPhiF)
{

/* For the ZPhiF doubles */

        *mem_ZPhiF=(ZPhiF *)realloc(*mem_ZPhiF, sofmZPhiF*sizeof(ZPhiF));
        if(!mem_ZPhiF )
        {
                printf( "failed to allocate memory for ZPhiF doubles\n ");
                exit(1);
        }
	memset(*mem_ZPhiF,0,sofmZPhiF*sizeof(ZPhiF));

  	return 1; 
}
