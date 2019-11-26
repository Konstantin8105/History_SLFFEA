/*
    This utility function allocates the memory for
    the truss finite element graphics program.

        Updated 3/16/00

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
#include "tsstrcgr.h"
#include "../truss/tsconst.h"
#include "../truss/tsstruct.h"

extern int input_flag, post_flag;

int tsMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, int sofmISTRESS )
{
/* Allocate the memory for when first running the program

        Updated 12/7/99

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

  	return 1; 
}

int tsMemory2_gr( XYZF **mem_XYZF, int sofmXYZF )
{
        *mem_XYZF=(XYZF *)calloc(sofmXYZF,sizeof(XYZF));
        if(!mem_XYZF )
        {
                printf( "failed to allocate memory for XYZF doubles\n ");
                exit(1);
        }
  	return 1; 
}

int tsReGetMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, int sofmISTRESS )
{
/* Allocate the memory for when first running the program

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

  	return 1; 
}

int tsReGetMemory2_gr( XYZF **mem_XYZF, int sofmXYZF )
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
