/*
    This program reads in the new input file and prepares it
    for graphical display.
  
   			Last Update 10/10/06

    SLFFEA source file
    Version:  1.1
    Copyright (C) 1999  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
 */

#if WINDOWS
#include <windows.h>
#endif

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../shell/shconst.h"
#include "../shell/shstruct.h"
#include "shstrcgr.h"
#include "../../common_gr/control.h"

/* glut header files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/********************* These are all the subroutines **************************/

/******** Data management and calculations ********/

void shforce_vectors0(int , BOUND , double *, XYZPhiF *);

void shdisp_vectors0(int , BOUND , double *);

void agvMakeAxesList(GLuint);

int shset( BOUND, int *, double *, double *, XYZPhiF *, SDIM *,
        ISTRAIN *, SDIM *, ISTRESS *, double *U, int * );

int shparameter( double *, SDIM *, SDIM *, double * );

int dotX(double *,double *, double *, int );

int shnormal_vectors (int *, double *, NORM * );

int shReGetMemory2_gr( XYZPhiF **, int );

int shreader_gr( FILE *, SDIM *, SDIM *);

int shreader( BOUND , int *, double *, int *, double *, MATL *,
        char *, FILE *, STRESS *, SDIM *, double *);

int shReGetMemory_gr( ISTRAIN **, ISTRESS **, NORM **, int, int );

int shReGetMemory( double **, int , int **, int , MATL **, int , XYZPhiI **,
	int , SDIM **, int, STRAIN **, STRESS **, int  );

int filecheck( char *, char *, FILE **, FILE **, FILE **, char * );

/******************************* GLOBAL VARIABLES **************************/

/****** FEA globals ******/
extern int dof, nmat, nmode, numel, numnp, integ_flag;
extern int stress_read_flag, element_stress_read_flag;
extern XYZPhiI *mem_XYZPhiI;
extern XYZPhiF *mem_XYZPhiF;
extern int *mem_int;
extern double *mem_double;
extern SDIM *mem_SDIM;
extern NORM *mem_NORM;
extern double *coord, *coord0;
extern double *U, *Uz_fib;
extern int *connecter;
extern BOUND bc;
extern MATL *matl;
extern int *el_matl;
extern double *force, *Voln;
extern STRESS *stress;
extern STRAIN *strain;
extern SDIM *stress_node;
extern SDIM *strain_node;

/* Global variables for the mesh color and nodal data */

extern ISTRESS *stress_color;
extern ISTRAIN *strain_color;
extern int *U_color, *el_matl_color;
extern NORM *norm, *norm0;
extern MATL *matl_crtl;

/* Global variables for drawing the axes */
extern GLuint AxesList, DispList, ForceList;   /* Display lists */

/* Global variables for drawing the force vectors */
extern XYZPhiF *force_vec, *force_vec0;

extern int input_flag, post_flag, matl_choice, node_choice, ele_choice, mode_choice;
extern int Before_flag, After_flag;
extern double amplify_factor;

int shGetNewMesh(void)
{
        int i, j, check;
        char *ccheck;
	int dum, dum1, dum2, dum3, dum4;
        double fpointx, fpointy, fpointz, node_vec[3], fdum;
	int sofmi, sofmf, sofmSTRESS, sofmISTRESS, sofmSTRAIN,
		sofmXYZPhiI, sofmXYZPhiF, sofmSDIM, sofmNORM, ptr_inc;
        char name[30], name2[30], osh_exten[4], buf[ BUFSIZ ];
        FILE *o1, *o2, *o3;

/* Delete the old display lists */

	glDeleteLists(AxesList,1);
	if(input_flag)
	{
		glDeleteLists(DispList,1);
		glDeleteLists(ForceList,1);
	}

/* Initialize filenames */

	memset(name,0,30*sizeof(char));
	memset(name2,0,30*sizeof(char));
	memset(osh_exten,0,4*sizeof(char));

/* Initialize old variables */

	input_flag = 1;
	post_flag = 1;
	After_flag = 1;
	Before_flag = 0;
	mode_choice = 0;
	amplify_factor = 1.0;
	matl_choice = 0;
	node_choice = 0;
	ele_choice = 0;

	ccheck = strncpy(osh_exten,".osh",4);
	if(!ccheck) printf( " Problems with strncpy \n");

        printf("What is the name of the input file containing the \n");
        printf("shell structural data? \n");
        scanf( "%30s",name2);

/*   o1 contains all the structural data for input
     o3 contains all the structural data for postprocessing
     o2 is used to determine the existance of input and post files
*/
        o2 = fopen( name2,"r" );
        if(o2 == NULL ) {
                printf("Can't find file %30s\n", name2);
                exit(1);
        }
        /*printf( "%3d %30s\n ",name2_length,name2);*/

        fgets( buf, BUFSIZ, o2 );
        fscanf( o2, "%d %d %d %d %d\n ",&numel,&numnp,&nmat,&nmode,&integ_flag);
        dof=numnp*ndof;
	nmode = abs(nmode);

/* Begin exmaining and checking for the existence of data files */

        check = filecheck( name, name2, &o1, &o2, &o3, osh_exten );
	if(!check) printf( " Problems with filecheck \n");

	if( input_flag )
	{
        	fgets( buf, BUFSIZ, o1 );
        	fscanf( o1, "%d %d %d %d %d\n ",&dum,&dum1,&dum2,&dum3,&dum4);
        	printf( "%d %d %d %d %d\n ",dum,dum1,dum2,dum3,dum4);
                /*printf( "name %30s\n ",name);*/
	}
	if( post_flag )
	{
        	fgets( buf, BUFSIZ, o3 );
        	fscanf( o3, "%d %d %d %d %d\n ",&dum,&dum1,&dum2,&dum3,&dum4);
        	printf( "%d %d %d %d %d\n ",dum,dum1,dum2,dum3,dum4);
                /*printf( "out %30s\n ",out);*/
	}

/*   Begin allocation of meomory */

/* For the doubles */
        sofmf=4*numnp*nsd+2*dof+numnp;

/* For the integers */
        sofmi= numel*npell+numel+numnp+1+1+dof;

/* For the XYZPhiI integers */
        sofmXYZPhiI=numnp+1+1;

/* For the SDIM doubles */
        sofmSDIM = 2*2*numnp;

/* For the STRESS */
	sofmSTRESS=1;

/* For the ISTRESS */
	sofmISTRESS=numel;

/* For the NORMS */
	sofmNORM=numel;
	if( input_flag && post_flag ) sofmNORM=2*numel;

        check = shReGetMemory( &mem_double, sofmf, &mem_int, sofmi, &matl, nmat,
                &mem_XYZPhiI, sofmXYZPhiI, &mem_SDIM, sofmSDIM, &strain,
		&stress, sofmSTRESS );
        if(!check) printf( " Problems with shReGetMemory \n");

        check = shReGetMemory_gr( &strain_color, &stress_color, &mem_NORM,
		sofmISTRESS, sofmNORM );
        if(!check) printf( " Problems with shReGetMemory_gr \n");

/* For the doubles */
                                        ptr_inc=0;
        coord=(mem_double+ptr_inc);     ptr_inc += 2*numnp*nsd;
        coord0=(mem_double+ptr_inc);    ptr_inc += 2*numnp*nsd;
        force=(mem_double+ptr_inc);     ptr_inc += dof;
        U=(mem_double+ptr_inc);         ptr_inc += dof;
        Uz_fib=(mem_double+ptr_inc);    ptr_inc += numnp;

/* For the materials */

	matl_crtl = matl;

/* For the integers */
                                                ptr_inc = 0;
        connecter=(mem_int+ptr_inc);            ptr_inc += numel*npell;
        el_matl=(mem_int+ptr_inc);              ptr_inc += numel;
        bc.force =(mem_int+ptr_inc);            ptr_inc += numnp+1;
        bc.num_force=(mem_int+ptr_inc);         ptr_inc += 1;
        U_color=(mem_int+ptr_inc);              ptr_inc += dof;

	el_matl_color = el_matl;

/* For the XYZPhiI integers */
                                          	ptr_inc = 0;
        bc.fix =(mem_XYZPhiI+ptr_inc);       	ptr_inc += numnp+1;
        bc.num_fix=(mem_XYZPhiI+ptr_inc);    	ptr_inc += 1;

/* For the SDIM doubles */
                                                ptr_inc = 0;
	stress_node=(mem_SDIM+ptr_inc);         ptr_inc += 2*numnp;
	strain_node=(mem_SDIM+ptr_inc);         ptr_inc += 2*numnp;

/* For the NORM doubles */
                                                ptr_inc = 0;
	norm =(mem_NORM+ptr_inc);
	if( input_flag && post_flag )           ptr_inc += numel;
	norm0 =(mem_NORM+ptr_inc);              ptr_inc += numel;

/* If there is no post file, then set coord, norm to coord0, norm0 */

	if( !post_flag )
	{
	    coord = coord0;
	    After_flag = 0;
	    Before_flag = 1;
	}

/* If there is no input file, then set coord0, norm0 to coord, and norm */

	if( !input_flag )
	{
	    /*coord0 = coord;*/
	    After_flag = 1;
	    Before_flag = 0;
	}

	stress_read_flag = 1;
	element_stress_read_flag = 0;
	if( post_flag )
	{
		check = shreader( bc, connecter, coord, el_matl, force, matl,
			name, o3, stress, stress_node, U);
        	if(!check) printf( " Problems with shreader \n");
		stress_read_flag = 0;

        	check = shreader_gr( o3, strain_node, stress_node);
        	if(!check) printf( " Problems with shreader_gr \n");
	}

	if( input_flag )
	{
		check = shreader( bc, connecter, coord0, el_matl, force, matl,
			name, o1, stress, stress_node, U);
        	if(!check) printf( " Problems with shreader \n");
	}

	if( post_flag )
	{
		check = shnormal_vectors(connecter, coord, norm );
                if(!check) printf( " Problems with shnormal_vectors \n");
	}

	if( input_flag )
	{
                check = shnormal_vectors(connecter, coord0, norm0 );
                if(!check) printf( " Problems with shnormal_vectors \n");
	}

	if( !input_flag )
	{
/* Calcuate the local z fiber displacement on each node */

            for( i = 0; i < numnp; ++i )
            {
/* Calcuate the projection of displacement vector to local z fiber vector */

	   	*(node_vec) = *(coord+nsd*(numnp+i))-*(coord+nsd*i);
	   	*(node_vec+1) = *(coord+nsd*(numnp+i)+1)-*(coord+nsd*i+1);
	   	*(node_vec+2) = *(coord+nsd*(numnp+i)+2)-*(coord+nsd*i+2);
	   	fdum = *(node_vec)*(*(node_vec)) +
		    *(node_vec+1)*(*(node_vec+1)) + *(node_vec+2)*(*(node_vec+2));
	   	fdum = sqrt(fdum);
	   	*(node_vec) /= fdum;
	   	*(node_vec+1) /= fdum;
	   	*(node_vec+2) /= fdum;
	   	check = check=dotX( (Uz_fib+i), (node_vec), (U+ndof*i), nsd);
	   	if(!check) printf( " Problems with dotX \n");

		*(coord0 + nsd*i) =
			*(coord+nsd*i) -  (*(U+ndof*i) +
			*(Uz_fib + i)*(*(U+ndof*i+4)));
		*(coord0 + nsd*i+1) =
			*(coord+nsd*i+1) -  (*(U+ndof*i+1) -
			*(Uz_fib + i)*(*(U+ndof*i+3)));
		*(coord0 + nsd*i+2) =
			*(coord+nsd*i+2) - *(U+ndof*i+2);

		*(coord0 + nsd*(numnp+i)) =
			*(coord+nsd*(numnp+i)) -  (*(U+ndof*i) +
			*(Uz_fib + i)*(*(U+ndof*i+4)));
		*(coord0 + nsd*(numnp+i) + 1) =
			*(coord+nsd*(numnp+i)+1) -  (*(U+ndof*i+1) -
			*(Uz_fib + i)*(*(U+ndof*i+3)));
		*(coord0 + nsd*(numnp+i) + 2) =
			*(coord+nsd*(numnp+i)+2) - *(U+ndof*i+2);
            }

	    /*for ( i = 0; i < numnp; ++i)
	    {
		*(coord0 + nsd*i) = *(coord+nsd*i) - *(U+ndof*i) -
			*(U+ndof*i+2)*(*(U+ndof*i+4));
		*(coord0 + nsd*i + 1) = *(coord+nsd*i+1) - *(U+ndof*i+1) +
			*(U+ndof*i+2)*(*(U+ndof*i+3));
		*(coord0 + nsd*i + 2) = *(coord+nsd*i+2) - *(U+ndof*i+2);
	    }*/
	}

/* Calcuate the local z fiber displacement on each node */

        for( i = 0; i < numnp; ++i )
        {
/* Calcuate the projection of displacement vector to local z fiber vector */

	   *(node_vec) = *(coord+nsd*(numnp+i))-*(coord+nsd*i);
	   *(node_vec+1) = *(coord+nsd*(numnp+i)+1)-*(coord+nsd*i+1);
	   *(node_vec+2) = *(coord+nsd*(numnp+i)+2)-*(coord+nsd*i+2);
	   fdum = *(node_vec)*(*(node_vec)) +
		*(node_vec+1)*(*(node_vec+1)) + *(node_vec+2)*(*(node_vec+2));
	   fdum = sqrt(fdum);
	   *(node_vec) /= fdum;
	   *(node_vec+1) /= fdum;
	   *(node_vec+2) /= fdum;
	   check = check=dotX( (Uz_fib+i), (node_vec), (U+ndof*i), nsd);
	   if(!check) printf( " Problems with dotX \n");

	   /*printf("\n node %3d %14.9f %14.9f %14.9f %14.9f %14.9f",
		i,*(Uz_fib+i),fdum,*(node_vec), *(node_vec+1), *(node_vec+2));*/
        }

/* For the XYZPhiF doubles */
        sofmXYZPhiF=2*bc.num_force[0];

        check = shReGetMemory2_gr( &mem_XYZPhiF, sofmXYZPhiF );
        if(!check) printf( " Problems with shReGetMemory2_gr \n");

                                                   ptr_inc = 0;
        force_vec =(mem_XYZPhiF+ptr_inc);          ptr_inc += bc.num_force[0];
        force_vec0 =(mem_XYZPhiF+ptr_inc);         ptr_inc += bc.num_force[0];

/* Search for extreme values */
 
/* In mesh viewer, search for extreme values of nodal points, displacements
   and stresss and strains to obtain viewing parameters and make color
   assignments.  Also initialize variables */

	check = shparameter( coord, strain_node, stress_node, U );
        if(!check) printf( " Problems with shparameter \n");

	check = shset( bc, connecter, coord, force, force_vec0, strain_node,
		strain_color, stress_node, stress_color, U, U_color );
        if(!check) printf( " Problems with shset \n");

  	AxesList = glGenLists(1);
  	agvMakeAxesList(AxesList);

	if( input_flag )
	{

/* create display list for displacement and force grapics vectors
   on undeformed mesh*/

  	    DispList = glGenLists(1);
  	    shdisp_vectors0(DispList, bc, coord0);

            for( i = 0; i < bc.num_force[0]; ++i)
            {
                fpointx = *(coord0+nsd*bc.force[i]);
                fpointy = *(coord0+nsd*bc.force[i] + 1);
                fpointz = *(coord0+nsd*bc.force[i] + 2);
                force_vec[i].x = fpointx - force_vec0[i].x;
                force_vec[i].y = fpointy - force_vec0[i].y;
                force_vec[i].z = fpointz - force_vec0[i].z;
                force_vec[i].phix = fpointx - force_vec0[i].phix;
		force_vec[i].phiy = fpointy - force_vec0[i].phiy;
            }
    
  	    ForceList = glGenLists(1);
  	    shforce_vectors0(ForceList, bc, coord0, force_vec);
    
	}

	if( post_flag )
	{
/* create force grapics vectors for deformed mesh*/

            for( i = 0; i < bc.num_force[0]; ++i)
            {
                fpointx = *(coord+nsd*bc.force[i]);
                fpointy = *(coord+nsd*bc.force[i] + 1);
                fpointz = *(coord+nsd*bc.force[i] + 2);
		force_vec[i].x = fpointx - force_vec0[i].x;
		force_vec[i].y = fpointy - force_vec0[i].y;
		force_vec[i].z = fpointz - force_vec0[i].z;
                force_vec[i].phix = fpointx - force_vec0[i].phix;
		force_vec[i].phiy = fpointy - force_vec0[i].phiy;
            }
	}

	return 1;
}
