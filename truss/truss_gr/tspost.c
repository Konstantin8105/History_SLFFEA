/*
    This program shows the 3 dimensional model of the finite
    element mesh for truss elements.
  
   			Last Update 12/4/06

    SLFFEA source file
    Version:  1.5
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006  San Le

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
 */

#if WINDOWS
#include <windows.h>
#endif

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "../truss/tsconst.h"
#include "../truss/tsstruct.h"
#include "tsgui.h"
#include "tsstrcgr.h"
#include "../../common_gr/control.h"
#include "../../common_gr/color_gr.h"

/* glut header files */
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#if LINUX
/* X11 header files. */
#include <GL/glx.h>

#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xmu/StdCmap.h>
#include <X11/keysym.h>

#include <X11/X.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>
#endif

/* The header files below can be used to incorporate
   widgets into your code */
#if 0
/* X11 Widget header files. */
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Form.h>
#include <GL/xmesa.h>
#include <GL/MesaDrawingArea.h>

#include <Xm/PushB.h>
#include <Xm/Form.h>
#include <GL/xmesa.h>
#include <GL/MesaMDrawingArea.h>

#define GLwMakeCurrent GLwMMakeCurrent
#endif

/********************* These are all the subroutines **************************/

/******** Data management and calculations ********/

void force_vectors0(int , BOUND , double *, XYZF *);

void disp_vectors0(int , BOUND , double *);

void agvMakeAxesList(GLuint);

int tsset( BOUND , double *, XYZF *, SDIM *, ISDIM *,
	SDIM *, ISDIM *, double *, int *);

int tsparameter( double *, SDIM *, SDIM *, double * );

int tsMemory2_gr( XYZF **, int );

int tsreader_gr( FILE *, SDIM *);

int tsreader( BOUND , int *, double *, int *, double *, MATL *,
	FILE *, SDIM *, double *);

int tsMemory_gr( ISDIM **, ISDIM **, int );

int tsMemory( double **, int, int **, int, MATL **, int , XYZI **, int,
	SDIM **, SDIM **, int );

int filecheck( char *, char *, FILE **, FILE **, FILE **, char *, int );

/******** For the Mesh ********/

void MeshKey_Special(int , int , int );

void tsMeshKeys( unsigned char , int , int  );

void tsmeshdraw(void);

void tsrender(void);

void MeshReshape(int , int );

void tsMeshDisplay(void);

void MeshInit(void);

/******** For the Control Panel ********/

void tsControlMouse(int , int , int , int );

void ControlReshape(int , int );

void tsControlDisplay(void);

void tsMenu();

void ControlInit(void);

int ControlDimInit( void );

/******** For the Screen Dump ********/

void ScreenShot(int, int);

/* These are functions provided by Phillip Winston for
   movement of the mesh.*/

void agvHandleButton(int , int , int , int );

void agvHandleMotion(int , int );


/******************************* GLOBAL VARIABLES **************************/

/****** FEA globals ******/
int dof, sdof, nmat, nmode, numel, numnp;
int stress_read_flag, element_stress_read_flag;
XYZI *mem_XYZI;
XYZF *mem_XYZF;
int *mem_int;
double *mem_double;
double *coord, *coord0;
double *U;
int *connecter;
BOUND bc;
MATL *matl;
int *el_matl;
double *force;
SDIM *stress;
SDIM *strain;

/* Global variables for the mesh color and nodal data */

ISDIM *stress_color;
ISDIM *strain_color;
int *U_color, *el_matl_color;
MATL *matl_crtl;

/* This used to be in MeshInit */

GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess[] = { 50.0 };
GLfloat light_position[] = { 1.0, 1.0, 6.0, 0.0 };

/****** graphics globals ******/

int choice_stress = 0;
int ControlWindow, MeshWindow;

/* Determines the step size for keyboard and Control panel movement */

double step_sizex = .1, step_sizey = .1, step_sizez = .1;

/****** Translation Variables ******/
double left_right, up_down, in_out, left_right0, up_down0, in_out0;

/****** Rotation Variables ******/
double xAngle = 0.0, yAngle = 0.0, zAngle = 0.0;

/****** Cross Section Plane Translation Variables ******/
double cross_sec_left_right, cross_sec_up_down, cross_sec_in_out,
	cross_sec_left_right0, cross_sec_up_down0, cross_sec_in_out0;

/* Global variables for drawing the axes */
GLuint AxesList, DispList, ForceList;   /* Display lists */
double AxisMax_x, AxisMax_y, AxisMax_z,
	AxisMin_x, AxisMin_y, AxisMin_z,
	IAxisMin_x, IAxisMin_y, IAxisMin_z;
double AxisLength_x, AxisLength_y, AxisLength_z, AxisLength_max, AxisPoint_step;

/* Global variables for drawing the force vectors */
XYZF *force_vec, *force_vec0;

/****** For drawing the Mesh Window ******/
double left, right, top, bottom, near, far, fscale, coord_rescale;
double ortho_left, ortho_right, ortho_top, ortho_bottom,
	ortho_left0, ortho_right0, ortho_top0, ortho_bottom0;
int mesh_width, mesh_height;
int com_mesh_width0 = mesh_width0;
int com_mesh_height0 = mesh_height0;

/* These Variables are for the Control Panel */

int sofi = sizeof(int);
int row_number = rowdim;
int ControlDiv_y[rowdim + 2], ControlDiv_x[rowdim + 2];
int control_width, control_height;
int com_control_width0 = control_width0;
int com_control_height0 =  control_height0;
int current_width, current_height;

double ratio = scaleFactor, ratio2 = 1.0;
double ratio_width = 1.0;
double ratio_height = 1.0;

int boxMove_x, boxMove_y, boxTextMove_x, textMove_x, textMove_y[rowdim+2];
int Color_flag[rowdim];

char RotateData[3][25] = { "    0.00", "    0.00", "    0.00" };
char MoveData[3][25] = { "    0.00", "    0.00", "    0.00" };
char AmplifyData[25] = { " 1.000e+00"};
char BoxData[2*boxnumber+2][25] = { "", "", "", "", "", "", "", "",
	"", "", "", "", "", "", "", "", "", "" };
char BoxText[25];

int del_height = 0;
int del_width = 0;

double matl_choicef = 0, node_choicef = 0, ele_choicef = 0;

int textDiv_xa = textDiv_xa0;
int textDiv_xb = textDiv_xb0;
int boxTextMove_x = boxTextMove_x0;

/* These Variables partition the stresses, strains, and displacements */

double Ux_div[boxnumber+1], Uy_div[boxnumber+1], Uz_div[boxnumber+1];
SDIM stress_div[boxnumber+1];
SDIM strain_div[boxnumber+1];
double init_right, init_left, init_top,
	init_bottom, init_near, init_far, true_far, dim_max;
SDIM del_stress, max_stress, min_stress;
SDIM del_strain, max_strain, min_strain;
double max_Ux, min_Ux, del_Ux, max_Uy, min_Uy, del_Uy,
	max_Uz, min_Uz, del_Uz, absolute_max_U, absolute_max_coord;

/* These are the flags */

int input_flag = 1,          /* Tells whether an input file exists or not */
    post_flag = 1,           /* Tells whether a post file exists or not   */
    color_choice = 1,        /* Set to desired color analysis(range 1 - 24) */
    choice = 0,              /* currently unused */
    matl_choice = -1,        /* Set to which material to view  */
    node_choice = -1,        /* Set to which node to view  */
    ele_choice = -1,         /* Set to which element to view  */
    mode_choice = 0;         /* Set to the desired modal post file */

int input_color_flag = 0;    /* Used with input_flag to determine how to draw input mesh */
int ortho_redraw_flag = 0;   /* calls MeshReshape(currently not used) */
int Dead_flag = 1,           /* Flag which is currently unused on and off  */
    Perspective_flag = 1,    /* Selects between orthographic and perspecive views */
    Render_flag = 0,         /* Selects between rendering or analyses */
    AppliedDisp_flag = 0,    /* Turns applied on and off */
    AppliedForce_flag = 0,   /* Turns applied force on and off */
    Material_flag = 0,       /* Turns material on and off */
    Node_flag = 0,           /* Turns Node ID on and off */
    Element_flag = 0,        /* Turns Element ID on and off */
    Axes_flag = 0,           /* Turns Axes on and off  */
    Transparent_flag = 0,    /* Turns Element Transparency on and off  */
    CrossSection_flag = 0;   /* Turns CrossSection_flag on and off  */
int Before_flag = 0,         /* Turns before mesh on and off */
    After_flag = 1,          /* Turns after mesh on and off */
    Both_flag = 0,           /* Turns both before and after on and off */
    Amplify_flag = 0;        /* Turns Amplification on and off */

double amplify_factor = 1.0; /* Amplifies deformed mesh for better viewing */
double amplify_step = 0.1;   /* Value of amplification icrement */
double amplify_step0 = 0.1;  /* Value of initial calculated amplification icrement */

int stress_flag = 0,       /* Tells whether stress viewing is on or off */
    strain_flag = 0,       /* Tells whether strain viewing is on and off */
    stress_strain = 0,     /* Used with above 2 flags to determine how to draw Control Panel */
    disp_flag = 0,         /* Tells whether displacement viewing is on and off */
    angle_flag = 0;        /* Tells whether angle viewing is on and off */


int main(int argc, char** argv)
{

	int i, j, check;
	char *ccheck;
	int dum, dum1, dum2, dum3;
	double fpointx, fpointy, fpointz;
	int sofmi, sofmf, sofmSDIM, sofmISDIM, sofmXYZI, sofmXYZF, ptr_inc;
	char name[30], name2[30], ots_exten[4], buf[ BUFSIZ ];
	int ots_exten_length = 4;
	FILE *o1, *o2, *o3;

	right=0;
	top=0;
	left=1000;
	bottom=1000;
	fscale=0;
	near=1.0;
	far=10.0;
	/*mesh_width=500;
	mesh_height=20;
	control_width=1000;
	control_height=1500;*/

/* Initialize filenames */

	memset(name,0,30*sizeof(char));
	memset(name2,0,30*sizeof(char));
	memset(ots_exten,0,ots_exten_length*sizeof(char));

	ccheck = strncpy(ots_exten,".ots",ots_exten_length);
	if(!ccheck) printf( " Problems with strncpy \n");

	printf("What is the name of the input file containing the \n");
	printf("truss structural data? (example: tower)\n");
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
	sdof=numnp*nsd;
	nmode = abs(nmode);

	check = filecheck( name, name2, &o1, &o2, &o3, ots_exten, ots_exten_length );
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
	sofmf=2*sdof+2*dof;

/* For the integers */
	sofmi= numel*npel+numel+numnp+1+1+dof;

/* For the XYZI integers */
	sofmXYZI=numnp+1+1;

/* For the SDIM */
	sofmSDIM=numel;

/* For the ISDIM */
	sofmISDIM=numel;

	check = tsMemory( &mem_double, sofmf, &mem_int, sofmi, &matl, nmat,
		&mem_XYZI, sofmXYZI, &strain, &stress, sofmSDIM );
	if(!check) printf( " Problems with tsMemory \n");

	check = tsMemory_gr( &strain_color, &stress_color, sofmISDIM );
	if(!check) printf( " Problems with tsMemory_gr \n");

/* For the doubles */
	                                ptr_inc=0;
	coord=(mem_double+ptr_inc);     ptr_inc += sdof;
	coord0=(mem_double+ptr_inc);    ptr_inc += sdof;
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

/* For the XYZI integers */
	                                        ptr_inc = 0;
	bc.fix =(mem_XYZI+ptr_inc);             ptr_inc += numnp+1;
	bc.num_fix=(mem_XYZI+ptr_inc);          ptr_inc += 1;


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
		check = tsreader( bc, connecter, coord, el_matl, force, matl,
			o3, stress, U);
		if(!check) printf( " Problems with tsreader \n");
		stress_read_flag = 0;

		check = tsreader_gr( o3, strain);
		if(!check) printf( " Problems with tsreader_gr \n");
	}
	if( input_flag )
	{
		check = tsreader( bc, connecter, coord0, el_matl, force, matl,
			o1, stress, U);
		if(!check) printf( " Problems with tsreader \n");

	}

/* For the XYZF doubles */
	sofmXYZF=2*bc.num_force[0];
/*
   This is allocated seperately from tsMemory_gr because we need to know the
   number of force vectors read from tsreader and stored in bc.num_force[0].
*/

	check = tsMemory2_gr( &mem_XYZF, sofmXYZF );
	if(!check) printf( " Problems with tsMemory2_gr \n");

	                                        ptr_inc = 0;
	force_vec =(mem_XYZF+ptr_inc);          ptr_inc += bc.num_force[0];
	force_vec0 =(mem_XYZF+ptr_inc);         ptr_inc += bc.num_force[0];

/* Search for extreme values */
 
/* In mesh viewer, search for extreme values of nodal points, displacements
   and stress and strains to obtain viewing parameters and make color
   assignments.  Also initialize variables */

	check = tsparameter( coord, strain, stress, U );
	if(!check) printf( " Problems with tsparameter \n");

/* Rescale undeformed coordinates */

	if( coord_rescale > 1.01 || coord_rescale < .99 )
	{
	   if( input_flag && post_flag )
	   {
		for( i = 0; i < numnp; ++i )
		{
			*(coord0+nsd*i) /= coord_rescale;
			*(coord0+nsd*i+1) /= coord_rescale;
			*(coord0+nsd*i+2) /= coord_rescale;
		}
	   }
	}

	if( !input_flag )
	{
	    for ( i = 0; i < numnp; ++i)
	    {
		*(coord0 + nsd*i) = *(coord+nsd*i) - *(U+ndof*i);
		*(coord0 + nsd*i + 1) = *(coord+nsd*i+1) - *(U+ndof*i+1);
		*(coord0 + nsd*i + 2) = *(coord+nsd*i+2) - *(U+ndof*i+2);
	    }
	}

	check = tsset( bc, force, force_vec0, strain, strain_color,
		stress, stress_color, U, U_color);
	if(!check) printf( " Problems with tsparameter \n");

/* Initialize the mesh viewer */

	glutInit(&argc, argv);
	glutInitWindowSize(mesh_width0, mesh_height0);
	glutInitWindowPosition(400, 215);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	MeshWindow = glutCreateWindow("SLFFEA");
	MeshInit ();

	AxesList = glGenLists(1);
	agvMakeAxesList(AxesList);

/* Below, I calculate force_vec[i].* for the force vectors graphics.  The reason I
   have coded things like this is because I think it gives me a slight improvement in
   speed.  When glutDisplayFunc displays the mesh, it continuously calls all the
   functions used in displaying the mesh like the subroutines which draw the force and
   prescribed displacement vectors.  This doesn't matter for the undeformed mesh where
   everything is drawn from display lists, but for the deformed mesh, it is an issue.
   So I calculate force_vec[i].* outside those functions, rather than simply passing
   force_vec0[i].* to the particular *force_vectors function and doing something like:

                fx = fpointx - force_vec0[node_num].x;
                fy = fpointy - force_vec0[node_num].y;
                fz = fpointz - force_vec0[node_num].z;

   There is probably only a small advantage, but that is the reason.
*/

	if( input_flag )
	{

/* create display list for displacement and force grapics vectors
   on undeformed mesh*/

	    DispList = glGenLists(1);
	    disp_vectors0(DispList, bc, coord0);

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
	    force_vectors0(ForceList, bc, coord0, force_vec);
    
	}

	if( post_flag )
	{
/* create force grapics vectors for undeformed mesh*/

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

/* Initiate variables in Control Panel */

	memset(Color_flag,0,rowdim*sofi);

	check = ControlDimInit();
	if(!check) printf( " Problems with ControlDimInit \n");

/* call display function  */

	glutDisplayFunc(tsMeshDisplay);

	glutReshapeFunc(MeshReshape);

/* Initialize Mouse Functions */

	glutMouseFunc(agvHandleButton);
	glutMotionFunc(agvHandleMotion);

/* Initialize Keyboard Functions */

	glutKeyboardFunc(tsMeshKeys);
	glutSpecialFunc(MeshKey_Special);

/* Initialize the Control Panel */

	glutInitWindowSize(control_width0, control_height0);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB );
	ControlWindow = glutCreateWindow("SLFFEA Control Panel");

	ControlInit();
	tsMenu();
	glutDisplayFunc(tsControlDisplay);
	glutReshapeFunc(ControlReshape);

	glutMouseFunc(tsControlMouse);

/* call function for hotkeys
 */
#if 0
	glutKeyboardFunc(ControlHandleKeys);
#endif
	glutMainLoop();

	free(strain);
	free(stress);
	free(strain_color);
	free(stress_color);
	free(matl);
	free(mem_double);
	free(mem_int);
	free(mem_XYZI);
	free(mem_XYZF);
	return 1;    /* ANSI C requires main to return int. */
}
