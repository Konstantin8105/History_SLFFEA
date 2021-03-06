/*
    This utility function allocates the memory for the graphics variables for
    the shell finite element graphics program.

        Updated 2/8/06

    SLFFEA source file
    Version:  1.5
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.

*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "shstrcgr.h"
#include "../shell/shconst.h"
#include "../shell/shstruct.h"

extern int input_flag, post_flag;

int shMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, int sofmISTRESS,
	NORM **mem_NORM, int sofmNORM )
{
/* Allocate the memory for when first running the program

	Updated 2/8/06

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

int shMemory2_gr( XYZPhiF **mem_XYZPhiF, int sofmXYZPhiF )
{
	*mem_XYZPhiF=(XYZPhiF *)calloc(sofmXYZPhiF,sizeof(XYZPhiF));
	if(!mem_XYZPhiF )
	{
		printf( "failed to allocate memory for XYZPhiF doubles\n ");
		exit(1);
	}
	return 1; 
}

int shReGetMemory_gr( ISTRAIN **strain_color, ISTRESS **stress_color, int sofmISTRESS,
	NORM **mem_NORM, int sofmNORM )
{
/* Allocate the memory for when first running the program

	Updated 2/8/06

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

int shReGetMemory2_gr( XYZPhiF **mem_XYZPhiF, int sofmXYZPhiF )
{
	*mem_XYZPhiF=(XYZPhiF *)realloc(*mem_XYZPhiF, sofmXYZPhiF*sizeof(XYZPhiF));
	if(!mem_XYZPhiF )
	{
		printf( "failed to allocate memory for XYZPhiF doubles\n ");
		exit(1);
	}
	memset(*mem_XYZPhiF,0,sofmXYZPhiF*sizeof(XYZPhiF));

	return 1; 
}
