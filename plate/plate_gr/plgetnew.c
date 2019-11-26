/* 
    This program reads in the new input file and prepares it
    for graphical display. 
  
   			Last Update 12/15/00

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
#include "../plate/plconst.h"
#include "../plate/plstruct.h"
#include "plstrcgr.h"
#include "../../common_gr/control.h"

/* glut header files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/********************* These are all the subroutines **************************/

/******** Data management and calculations ********/

void plforce_vectors0(int , BOUND , double *, ZPhiF *);

void pldisp_vectors0(int , BOUND , double *);

void agvMakeAxesList(GLuint);

int plset( BOUND bc, int *, double *, MDIM *, ICURVATURE *, double *, ZPhiF *,
        MDIM *, IMOMENT *, SDIM *, ISTRAIN *, SDIM *, ISTRESS *,
        double *, int * );

int plparameter( double *, MDIM *, MDIM *, SDIM *, SDIM *, double * );

int plnormal_vectors (int *, double *, NORM *, double * );

int plReGetMemory2_gr( ZPhiF **, int );

int plreader_gr( FILE *, MDIM *, MDIM *, SDIM *, SDIM *);

int plreader( BOUND , int *, double *, int *, double *, MATL *, MOMENT *,
	MDIM *, char *, FILE *, STRESS *, SDIM *, double *);

int plReGetMemory_gr( ICURVATURE **, IMOMENT **, ISTRAIN **, ISTRESS **, NORM **,
        int, int );

int plReGetMemory( double **, int , int **, int , MATL **, int , ZPhiI **, int ,
        MDIM **, SDIM **, int , CURVATURE **, MOMENT **, STRAIN **, STRESS **, int );

int filecheck( char *, char *, FILE **, FILE **, FILE **, char * );

/******************************* GLOBAL VARIABLES **************************/

/****** FEA globals ******/
extern int dof, nmat, nmode, numel, numnp;
extern int stress_read_flag, element_stress_read_flag;
extern ZPhiI *mem_ZPhiI;
extern ZPhiF *mem_ZPhiF;
extern int *mem_int;
extern double *mem_double;
extern SDIM *mem_SDIM;
extern MDIM *mem_MDIM;
extern NORM *mem_NORM;
extern double *coord, *coord0, *zcoord;
extern double *U;
extern int *connecter;
extern BOUND bc;
extern MATL *matl;
extern int *el_matl;
extern double *force;
extern MOMENT *moment;
extern STRESS *stress;
extern CURVATURE *curve;
extern STRAIN *strain;
extern MDIM *moment_node;
extern SDIM *stress_node;
extern MDIM *curve_node;
extern SDIM *strain_node;

/* Global variables for the mesh color and nodal data */

extern IMOMENT *moment_color;
extern ICURVATURE *curve_color;
extern ISTRESS *stress_color;
extern ISTRAIN *strain_color;
extern int *U_color, *el_matl_color;
extern NORM *norm, *norm0;
extern MATL *matl_crtl;

/* Global variables for drawing the axes */
extern GLuint AxesList, DispList, ForceList;   /* Display lists */

/* Global variables for drawing the force vectors */
extern ZPhiF *force_vec, *force_vec0;

extern int input_flag, post_flag, matl_choice, node_choice, ele_choice, mode_choice; 
extern int Before_flag, After_flag;
extern double amplify_factor;

int plGetNewMesh()
{
        int i, j, check;
        char *ccheck;
	int dum, dum1, dum2, dum3;
        double fpointx, fpointy, fpointz;
	int  sofmi, sofmf, sofmSTRESS, sofmISTRESS, sofmSTRAIN,
		sofmSDIM, sofmZPhiI, sofmZPhiF, sofmNORM, ptr_inc;
        char name[30], name2[30], opl_exten[4], buf[ BUFSIZ ];
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
	memset(opl_exten,0,4*sizeof(char));

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

	ccheck = strncpy(opl_exten,".opl",4);
	if(!ccheck) printf( " Problems with strncpy \n");

        printf("What is the name of the input file containing the \n");
        printf("plate structural data? \n");
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
        fscanf( o2, "%d %d %d %d\n ",&numel,&numnp,&nmat,&nmode);
        dof=numnp*ndof;
	nmode = abs(nmode);

/* Begin exmaining and checking for the existence of data files */

	check = filecheck( name, name2, &o1, &o2, &o3, opl_exten );
	if(!check) printf( " Problems with filecheck \n");

	if( input_flag )
	{
        	fgets( buf, BUFSIZ, o1 );
        	fscanf( o1, "%d %d %d %d\n ",&dum,&dum1,&dum2,&dum3);
        	printf( "%d %d %d %d\n ",dum,dum1,dum2,dum3);
                /*printf( "name %30s\n ",name);*/
	}
	if( post_flag )
	{
        	fgets( buf, BUFSIZ, o3 );
        	fscanf( o3, "%d %d %d %d\n ",&dum,&dum1,&dum2,&dum3);
        	printf( "%d %d %d %d\n ",dum,dum1,dum2,dum3);
                /*printf( "out %30s\n ",out);*/
	}

/*   Begin allocation of meomory */

/* For the doubles */
        sofmf=2*numnp*nsd+numnp+2*dof;

/* For the integers */
        sofmi= numel*npel+numel+numnp+1+1+dof;

/* For the ZPhiI integers */
        sofmZPhiI=numnp+1+1;

/* For the SDIM doubles */
        sofmSDIM = 2*numnp;

/* For the STRESS */
        sofmSTRESS=1;

/* For the ISTRESS */
        sofmISTRESS=numel;

/* For the NORMS */
	sofmNORM=numel;
	if( input_flag && post_flag ) sofmNORM=2*numel;

	check = plReGetMemory( &mem_double, sofmf, &mem_int, sofmi, &matl, nmat, &mem_ZPhiI,
		sofmZPhiI, &mem_MDIM, &mem_SDIM, sofmSDIM, &curve, &moment, &strain,
		&stress, sofmSTRESS );
	if(!check) printf( " Problems with plReGetMemory \n");

	check = plReGetMemory_gr( &curve_color, &moment_color, &strain_color,
		&stress_color, &mem_NORM, sofmISTRESS, sofmNORM );
	if(!check) printf( " Problems with plReGetMemory_gr \n");

/* For the doubles */
                                        ptr_inc=0;
        coord=(mem_double+ptr_inc);     ptr_inc += numnp*nsd;
        coord0=(mem_double+ptr_inc);    ptr_inc += numnp*nsd;
        zcoord=(mem_double+ptr_inc);    ptr_inc += numnp;
        force=(mem_double+ptr_inc);     ptr_inc += dof;
        U=(mem_double+ptr_inc);         ptr_inc += dof;

/* For the materials */

	matl_crtl = matl;

/* For the integers */
                                                ptr_inc = 0;
        connecter=(mem_int+ptr_inc);            ptr_inc += numel*npel;
        el_matl=(mem_int+ptr_inc);              ptr_inc += numel;
        bc.force =(mem_int+ptr_inc);            ptr_inc += numnp+1;
        bc.num_force=(mem_int+ptr_inc);         ptr_inc += 1;
        U_color=(mem_int+ptr_inc);              ptr_inc += dof;

	el_matl_color = el_matl;

/* For the ZPhiI integers */
                                          	ptr_inc = 0;
        bc.fix =(mem_ZPhiI+ptr_inc);       	ptr_inc += numnp+1;
        bc.num_fix=(mem_ZPhiI+ptr_inc);    	ptr_inc += 1;

/* For the SDIM doubles */
                                                ptr_inc = 0;
	stress_node=(mem_SDIM+ptr_inc);         ptr_inc += numnp;
	strain_node=(mem_SDIM+ptr_inc);         ptr_inc += numnp;

/* For the MDIM doubles */
                                                ptr_inc = 0;
	moment_node=(mem_MDIM+ptr_inc);         ptr_inc += numnp;
	curve_node=(mem_MDIM+ptr_inc);          ptr_inc += numnp;

/* For the NORM doubles */
                                                ptr_inc = 0;
	norm =(mem_NORM+ptr_inc);
	if( input_flag && post_flag )           ptr_inc += numel;
	norm0 =(mem_NORM+ptr_inc);              ptr_inc += numel;

/* If there is no post file, then set coord to coord0 */

	if( !post_flag )
	{
	    coord = coord0;
	    After_flag = 0;
	    Before_flag = 1;
	}

/* If there is no input file, then set coord0 to coord */

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
		check = plreader( bc, connecter, coord, el_matl, force, matl, moment,
			moment_node, name, o3, stress, stress_node, U);
        	if(!check) printf( " Problems with plreader \n");
		stress_read_flag = 0;

        	check = plreader_gr( o3, curve_node, moment_node, strain_node,
			stress_node);
        	if(!check) printf( " Problems with plreader_gr \n");

/* Set zcoord to U */

		for ( i = 0; i < numnp; ++i)
		{
			*(zcoord + i) = *(U+ndof*i);
		}
	}

	if( input_flag )
	{
		check = plreader( bc, connecter, coord0, el_matl, force, matl, moment,
			moment_node, name, o1, stress, stress_node, U);
        	if(!check) printf( " Problems with plreader \n");

/* Set zcoord to 0.0 */

		for ( i = 0; i < numnp; ++i)
		{
			*(zcoord + i) = 0.0;
		}
	}

	if( post_flag )
	{
		check = plnormal_vectors(connecter, coord, norm, zcoord );
                if(!check) printf( " Problems with plnormal_vectors \n");
	}

	if( input_flag )
	{
                check = plnormal_vectors(connecter, coord0, norm0, zcoord );
                if(!check) printf( " Problems with plnormal_vectors \n");
	}

	if( !input_flag )
	{
	    for ( i = 0; i < numnp; ++i)
	    {
		*(coord0 + nsd*i) = *(coord+nsd*i) - *(U+ndof*i)*(*(U+ndof*i+2));
		*(coord0 + nsd*i + 1) = *(coord+nsd*i+1) + *(U+ndof*i)*(*(U+ndof*i+1));

	    }
	}

/* Initialize zcoord */

	for ( i = 0; i < numnp; ++i)
	{
		*(zcoord + i) = *(U+ndof*i);
	}

/* For the ZPhiF doubles */
        sofmZPhiF=2*bc.num_force[0];

	check = plReGetMemory2_gr( &mem_ZPhiF, sofmZPhiF );
	if(!check) printf( " Problems with plReGetMemory2_gr \n");

                                                 ptr_inc = 0;
        force_vec =(mem_ZPhiF+ptr_inc);          ptr_inc += bc.num_force[0];
        force_vec0 =(mem_ZPhiF+ptr_inc);         ptr_inc += bc.num_force[0];

/* Search for extreme values */
 
/* In mesh viewer, search for extreme values of nodal points, displacements
   and stresss and strains to obtain viewing parameters and make color
   assignments.  Also initialize variables */

	check = plparameter( coord, curve_node, moment_node, strain_node,
		stress_node, U );
        if(!check) printf( " Problems with plparameter \n");

	check = plset( bc, connecter, coord, curve_node, curve_color, force,
		force_vec0, moment_node, moment_color, strain_node, strain_color,
		stress_node, stress_color, U, U_color );
        if(!check) printf( " Problems with plset \n");

  	AxesList = glGenLists(1);
  	agvMakeAxesList(AxesList);

	if( input_flag )
	{

/* create display list for displacement and force grapics vectors
   on undeformed mesh*/

  	    DispList = glGenLists(1);
  	    pldisp_vectors0(DispList, bc, coord0);

            for( i = 0; i < bc.num_force[0]; ++i)
            {
                fpointx = *(coord0+nsd*bc.force[i]);
                fpointy = *(coord0+nsd*bc.force[i] + 1);
                fpointz = 0.0;
                force_vec[i].z = fpointz - force_vec0[i].z;
                force_vec[i].phix = fpointx - force_vec0[i].phix;
		force_vec[i].phiy = fpointy - force_vec0[i].phiy;
            }
    
  	    ForceList = glGenLists(1);
  	    plforce_vectors0(ForceList, bc, coord0, force_vec);
    
	}

	if( post_flag )
	{
/* create force grapics vectors for deformed mesh*/

            for( i = 0; i < bc.num_force[0]; ++i)
            {
                fpointx = *(coord+nsd*bc.force[i]);
                fpointy = *(coord+nsd*bc.force[i] + 1);
                fpointz = *(zcoord+bc.force[i]);
		force_vec[i].z = fpointz - force_vec0[i].z;
                force_vec[i].phix = fpointx - force_vec0[i].phix;
		force_vec[i].phiy = fpointy - force_vec0[i].phiy;
            }
	}

	return 1;
}
