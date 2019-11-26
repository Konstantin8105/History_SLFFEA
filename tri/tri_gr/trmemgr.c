/*
    This utility function allocates the memory for
    the triangle finite element graphics program.

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
#include "trstrcgr.h"
#include "../tri/trconst.h"
#include "../tri/trstruct.h"

extern int input_flag, post_flag;

int trMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, int sofmISTRESS )
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

int trMemory2_gr( XYF **mem_XYF, int sofmXYF )
{
	*mem_XYF=(XYF *)calloc(sofmXYF,sizeof(XYF));
	if(!mem_XYF )
	{
		printf( "failed to allocate memory for XYF doubles\n ");
		exit(1);
	}
	return 1;
}

int trReGetMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, int sofmISTRESS )
{
/* Re-Allocate the memory for new input file

        Updated 12/23/99

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

int trReGetMemory2_gr( XYF **mem_XYF, int sofmXYF )
{
	*mem_XYF=(XYF *)realloc(*mem_XYF, sofmXYF*sizeof(XYF));
	if(!mem_XYF )
	{
		printf( "failed to allocate memory for XYF doubles\n ");
		exit(1);
	}
	memset(*mem_XYF,0,sofmXYF*sizeof(XYF));

	return 1;
}

