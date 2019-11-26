/*
    This utility function allocates the memory for
    the brick finite element graphics program. 

        Updated 5/14/00

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.

*/

#include <stdlib.h>
#include <string.h>
#include "brstrcgr.h"
#include "../brick/brconst.h"
#if BRICK1
#include "../brick/brstruct.h"
#endif
#if BRICK2
#include "../brick2/br2struct.h"
#endif

extern int input_flag, post_flag;

int brMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, NORM **mem_NORM,
	int sofmISTRESS, int sofmNORM )
{
/* Allocate the memory for when first running the program

        Updated 12/23/99

*/

/* For the ISTRESS integers */

        *stress_color=(ISTRESS *)calloc(sofmISTRESS,sizeof(ISTRESS));
        if(!stress_color )
        {
                printf( "failed to allocate memory for stress integers\n ");
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

int brMemory2_gr( XYZF **mem_XYZF, int sofmXYZF )
{
        *mem_XYZF=(XYZF *)calloc(sofmXYZF,sizeof(XYZF));
        if(!mem_XYZF )
        {
                printf( "failed to allocate memory for XYZF doubles\n ");
                exit(1);
        }
  	return 1; 
}


int brReGetMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, NORM **mem_NORM,
	int sofmISTRESS, int sofmNORM )
{
/* Re-Allocate the memory for new input file 

        Updated 3/14/00

*/

/* For the ISTRESS integers */

        *stress_color=(ISTRESS *)realloc(*stress_color, sofmISTRESS*sizeof(ISTRESS));
        if(!stress_color )
        {
                printf( "failed to allocate memory for stress integers\n ");
                exit(1);
        }
	memset(*stress_color,0,sofmISTRESS*sizeof(ISTRESS));

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

int brReGetMemory2_gr( XYZF **mem_XYZF, int sofmXYZF )
{
        *mem_XYZF=(XYZF *)realloc(*mem_XYZF, sofmXYZF*sizeof(XYZF));
        if(!mem_XYZF )
        {
                printf( "failed to allocate memory for XYZF doubles\n ");
                exit(1);
        }
	memset(*mem_XYZF,0,sofmXYZF*sizeof(XYZF));

  	return 1; 
}