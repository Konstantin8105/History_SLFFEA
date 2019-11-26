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
#include "../brick/brconst.h"
#include "../brick2/br2struct.h"
#include "../brick_gr/brstrcgr.h"
#include "../../common_gr/control.h"

/* glut header files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

/********************* These are all the subroutines **************************/

/******** Data management and calculations ********/

void brforce_vectors0(int , BOUND , double *, XYZF *);

void brdisp_vectors0(int , BOUND , double *);

void agvMakeAxesList(GLuint);

int br2set(BOUND , double *, int *, double *, XYZF *,  double *, int *, SDIM *,
	ISTRAIN *, SDIM *, ISTRESS *, double *, int *, double *, int *);

int br2parameter( double *, double *, SDIM *, SDIM *, double *, double *);

int brnormal_vectors (int *, double *, NORM * );

int brReGetMemory2_gr( XYZF **, int );

int brConnectSurfreader( int *, int *, char *);

int brreader_gr( FILE *, SDIM *, SDIM *);

int br2reader( BOUND , int *, int *, double *, int *, int *, double *, double *,
	double *, MATL *, char *name, FILE *, double *, STRESS *, SDIM *, double *,
	double *, double *);

int brReGetMemory_gr( ISTRAIN **, ISTRESS **, NORM **, int, int );

int brReGetMemory( double **, int, int **, int, MATL **, int , XYZI **, int,
	SDIM **, int, STRAIN **, STRESS **, int );

int filecheck( char *, char *, FILE **, FILE **, FILE **, char * );

/******************************* GLOBAL VARIABLES **************************/

/****** FEA globals ******/
extern int dof, Tdof, nmat, nmode, numel, numel_film, numnp;
extern int stress_read_flag, element_stress_read_flag;
extern int linear_analysis, thermal_analysis;
extern XYZI *mem_XYZI;
extern XYZF *mem_XYZF;
extern int *mem_int;
extern double *mem_double;
extern SDIM *mem_SDIM;
extern NORM *mem_NORM;
extern double *coord, *coord0;
extern double *U;
extern int *connecter, *connecter_film;
extern BOUND bc;
extern MATL *matl;
extern int *el_matl;
extern int *el_matl_film;
extern double *force;
extern STRESS *stress;
extern STRAIN *strain;
extern SDIM *stress_node;
extern SDIM *strain_node;

extern double *heat_el, *heat_node, *T, *TB, *Q;


/* Global variables for the mesh color and nodal data */

extern ISTRESS *stress_color;
extern ISTRAIN *strain_color;
extern int *U_color, *T_color, *Q_color, *el_matl_color;
extern NORM *norm, *norm0;
extern MATL *matl_crtl;

/* Global variables for drawing the axes */
extern GLuint AxesList, DispList, ForceList;   /* Display lists */

/* Global variables for drawing the force vectors */
extern XYZF *force_vec, *force_vec0;

extern int input_flag, post_flag, matl_choice, node_choice, ele_choice, mode_choice;
extern int Before_flag, After_flag;
extern double amplify_factor;

int br2GetNewMesh(void)
{
        int i, j, check;
        char *ccheck;
        int dum, dum1, dum2, dum3, dum4, dum5;
        double fpointx, fpointy, fpointz;
	int  sofmi, sofmf, sofmSTRESS, sofmISTRESS, sofmSTRAIN,
		sofmXYZI, sofmXYZF, sofmSDIM, sofmNORM, ptr_inc;
        char name[30], name2[30], obr_exten[4], buf[ BUFSIZ ];
        FILE *o1, *o2, *o3, *o4;

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
	memset(obr_exten,0,4*sizeof(char));

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

	ccheck = strncpy(obr_exten,".obr",4);
	if(!ccheck) printf( " Problems with strncpy \n");

#if LINUX
        printf("What is the name of the input file containing the \n");
        printf("brick 2 structural data? \n");
        scanf( "%30s",name2);
#endif

#if WINDOWS
        o4 = fopen( "brname.txt", "r");

        if(o4 == NULL ) {
                printf("Can't find file brname.txt\n");
                exit(1);
        }
        fscanf(o4, "%30s\n", name2);
	fclose(o4);
#endif

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
        fscanf( o2, "%d %d %d %d %d %d\n ",
		&numel,&numnp,&nmat,&numel_film,&linear_analysis,&thermal_analysis);
	Tdof=numnp*Tndof;
        dof=numnp*ndof;

/* Begin exmaining and checking for the existence of data files */

	check = filecheck( name, name2, &o1, &o2, &o3, obr_exten );
	if(!check) printf( " Problems with filecheck \n");

	if( input_flag )
	{
        	fgets( buf, BUFSIZ, o1 );
        	fscanf( o1, "%d %d %d %d %d %d\n ",&dum,&dum1,&dum2,&dum3,&dum4,&dum5);
        	printf( "%d %d %d %d %d %d\n ",dum,dum1,dum2,dum3,dum4,dum5);
                /*printf( "name %30s\n ",name);*/
	}
	if( post_flag )
	{
        	fgets( buf, BUFSIZ, o3 );
        	fscanf( o3, "%d %d %d %d %d %d\n ",&dum,&dum1,&dum2,&dum3,&dum4,&dum5);
        	printf( "%d %d %d %d %d %d\n ",dum,dum1,dum2,dum3,dum4,dum5);
                /*printf( "out %30s\n ",out);*/
	}

/*   Begin allocation of meomory */

/* For the doubles */
	sofmf=2*numnp*nsd + dof + numel + 4*Tdof + dof;

/* For the integers */
	sofmi= numel*npel + numel_film*npel_film + numel + numel_film +
		numnp + 1 + numel + 1 + 4*numnp + 4 + 6 + dof + 2*Tdof;

/* For the XYZI integers */
        sofmXYZI=numnp+1+1;

/* For the SDIM doubles */
	sofmSDIM = 2*numnp;

/* For the STRESS */
	sofmSTRESS=1;

/* For the ISTRESS */
	sofmISTRESS=numel;

/* For the NORMS */
	sofmNORM=numel;
	if( input_flag && post_flag ) sofmNORM=2*numel;

	check = brReGetMemory( &mem_double, sofmf, &mem_int, sofmi, &matl, nmat,
		&mem_XYZI, sofmXYZI, &mem_SDIM, sofmSDIM, &strain, &stress, sofmSTRESS );
	if(!check) printf( " Problems with brReGetMemory \n");

	check = brReGetMemory_gr( &strain_color, &stress_color, &mem_NORM, sofmISTRESS,
		sofmNORM );
	if(!check) printf( " Problems with brReGetMemory_gr \n");

/* For the doubles */
                                                ptr_inc=0;
	coord=(mem_double+ptr_inc);             ptr_inc += numnp*nsd;
	coord0=(mem_double+ptr_inc);            ptr_inc += numnp*nsd;
	force=(mem_double+ptr_inc);             ptr_inc += dof;
	heat_el=(mem_double+ptr_inc);           ptr_inc += numel;
	heat_node=(mem_double+ptr_inc);         ptr_inc += Tdof;
	Q=(mem_double+ptr_inc);                 ptr_inc += Tdof;
	T=(mem_double+ptr_inc);                 ptr_inc += Tdof;
	TB=(mem_double+ptr_inc);                ptr_inc += Tdof;
	U=(mem_double+ptr_inc);                 ptr_inc += dof;

/* For the materials */

	matl_crtl = matl;

/* For the integers */
                                                ptr_inc = 0;
	connecter=(mem_int+ptr_inc);            ptr_inc += numel*npel;
	connecter_film=(mem_int+ptr_inc);       ptr_inc += numel_film*npel_film;
	el_matl=(mem_int+ptr_inc);              ptr_inc += numel;
	el_matl_film=(mem_int+ptr_inc);         ptr_inc += numel_film;
	bc.force =(mem_int+ptr_inc);            ptr_inc += numnp+1;
	bc.heat_el =(mem_int+ptr_inc);          ptr_inc += numel+1;
	bc.heat_node =(mem_int+ptr_inc);        ptr_inc += numnp+1;
	bc.Q =(mem_int+ptr_inc);                ptr_inc += numnp+1;
	bc.T =(mem_int+ptr_inc);                ptr_inc += numnp+1;
	bc.TB =(mem_int+ptr_inc);               ptr_inc += numnp+1;
	bc.num_force=(mem_int+ptr_inc);         ptr_inc += 1;
	bc.num_heat_el=(mem_int+ptr_inc);       ptr_inc += 1;
	bc.num_heat_node=(mem_int+ptr_inc);     ptr_inc += 1;
	bc.num_Q=(mem_int+ptr_inc);             ptr_inc += 1;
	bc.num_T=(mem_int+ptr_inc);             ptr_inc += 1;
	bc.num_TB=(mem_int+ptr_inc);            ptr_inc += 1;
	U_color=(mem_int+ptr_inc);              ptr_inc += dof;
	T_color=(mem_int+ptr_inc);              ptr_inc += Tdof;
	Q_color=(mem_int+ptr_inc);              ptr_inc += Tdof;

	el_matl_color = el_matl;

/* For the XYZI integers */
                                          	ptr_inc = 0;
        bc.fix =(mem_XYZI+ptr_inc);       	ptr_inc += numnp+1;
        bc.num_fix=(mem_XYZI+ptr_inc);    	ptr_inc += 1;

/* For the SDIM doubles */
                                                ptr_inc = 0;
	stress_node=(mem_SDIM+ptr_inc);         ptr_inc += numnp;
	strain_node=(mem_SDIM+ptr_inc);         ptr_inc += numnp;

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
		check = br2reader( bc, connecter, connecter_film, coord, el_matl,
			el_matl_film, force, heat_el, heat_node, matl, name, o3,
			Q, stress, stress_node, T, TB, U);
        	if(!check) printf( " Problems with br2reader \n");
		stress_read_flag = 0;

        	check = brreader_gr( o3, strain_node, stress_node);
        	if(!check) printf( " Problems with brreader_gr \n");
	}
	if( input_flag )
	{
		check = br2reader( bc, connecter, connecter_film, coord0, el_matl,
			el_matl_film, force, heat_el, heat_node, matl, name, o1,
			Q, stress, stress_node, T, TB, U);
        	if(!check) printf( " Problems with br2reader \n");
	}

	check = brConnectSurfreader( connecter, el_matl, name);
	if(!check) printf( " Problems with brConnectSurfreader \n");

	if( post_flag )
	{
		check = brnormal_vectors(connecter, coord, norm );
        	if(!check) printf( " Problems with brnormal_vectors \n");
	}

	if( input_flag )
	{
		check = brnormal_vectors(connecter, coord0, norm0 );
        	if(!check) printf( " Problems with brnormal_vectors \n");
	}

	if( !input_flag )
	{
	    for ( i = 0; i < numnp; ++i)
	    {
	    	*(coord0 + nsd*i) = *(coord + nsd*i) - *(U + ndof*i);
	    	*(coord0 + nsd*i + 1) = *(coord + nsd*i + 1) - *(U + ndof*i + 1);
	    	*(coord0 + nsd*i + 2) = *(coord + nsd*i + 2) - *(U + ndof*i + 2);
	    }
	}

/* For the XYZF doubles */
        sofmXYZF=2*bc.num_force[0];

	check = brReGetMemory2_gr( &mem_XYZF, sofmXYZF );
	if(!check) printf( " Problems with brReGetMemory2_gr \n");

                                                ptr_inc = 0;
        force_vec =(mem_XYZF+ptr_inc);          ptr_inc += bc.num_force[0];
        force_vec0 =(mem_XYZF+ptr_inc);         ptr_inc += bc.num_force[0];

/* Search for extreme values */
 
/* In mesh viewer, search for extreme values of nodal points, displacements
   and stresss and strains to obtain viewing parameters and set color
   assignments.  Also initialize variables */

	check = br2parameter( coord, Q, strain_node, stress_node, T, U);
        if(!check) printf( " Problems with br2parameter \n");

	check = br2set( bc, coord, connecter, force, force_vec0, Q, Q_color, strain_node,
		strain_color, stress_node, stress_color, T, T_color, U, U_color);
        if(!check) printf( " Problems with br2set \n");

  	AxesList = glGenLists(1);
  	agvMakeAxesList(AxesList);

	if( input_flag )
	{

/* create display list for displacement and force grapics vectors
   on undeformed mesh*/

  	    DispList = glGenLists(1);
  	    brdisp_vectors0(DispList, bc, coord0);

            for( i = 0; i < bc.num_force[0]; ++i)
            {
                fpointx = *(coord0+nsd*bc.force[i]);
                fpointy = *(coord0+nsd*bc.force[i] + 1);
                fpointz = *(coord0+nsd*bc.force[i] + 2);
                force_vec[i].x = fpointx - force_vec0[i].x;
		force_vec[i].y = fpointy - force_vec0[i].y;
		force_vec[i].z = fpointz - force_vec0[i].z;
            }
    
  	    ForceList = glGenLists(1);
  	    brforce_vectors0(ForceList, bc, coord0, force_vec);
    
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
            }
	}

	return 1;
}
