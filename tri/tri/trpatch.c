/*
    This library function reads in data from a finite element
    data set to prepare it for the patch test analysis on a
    triangle element.  First, it reads in the data, then
    it creates the prescribed displacements so that the main
    finite element program can do the analysis. 

    There are 3 output patch test files: patch.x, patch.y, and patch.z.
    This triangle element is like the truss element in that it can not
    handle loads perpendicular to the local plane of an element because
    it is simply a membrane.  I had tried to do a mesh which had x, y, and
    z components with displacements colinear with the plane, but the
    results showed that this was not a valid analysis.  You can see this
    by taking a single element in the global x-y plane and leaving one node
    completely free.  Even if you put a load with only x-y components, the
    global stiffness will still be singular because the z DOF needs to be
    fixed.  If the calculation does not blow up, it will still produce
    incorrect results.  This is why the patch test needs to be broken down
    into 3 meshes with only 2 DOFs on the free center node corresponding to
    the plane of the mesh for each patch.
    

		Updated 10/25/05

    SLFFEA source file
    Version:  1.5
    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "trconst.h"
#include "trstruct.h"


int main(int argc, char** argv)
{
	int dof, nmat, numel, numnp, nmode, plane_stress_flag;
        int i,j,dum,dum2,dum3,dum4;
	double fdum1, fdum2, fdum3, ux2[9], uy2[9], uz2[9],
		ux3[9], uy3[9], uz3[9], ux4[9], uy4[9], uz4[9];
        char name[20];
	char buf[ BUFSIZ ];
        FILE *o1, *o2, *o3, *o4;
	char text;

        printf( "What is the name of the file containing the \n");
        printf( "triangle structural data? \n");
        scanf( "%20s",name);

        o1 = fopen( name,"r" );
        if(o1 == NULL ) {
                printf("error on open\n");
                exit(1);
        }
        o2 = fopen( "patch.x","w" );
        o3 = fopen( "patch.y","w" );
        o4 = fopen( "patch.z","w" );

        fprintf( o2, "   numel numnp nmat plane_stress_flag  This is the patch test \n ");
        fprintf( o3, "   numel numnp nmat plane_stress_flag  This is the patch test \n ");
        fprintf( o4, "   numel numnp nmat plane_stress_flag  This is the patch test \n ");
        fgets( buf, BUFSIZ, o1 );
        fscanf( o1, "%d %d %d %d %d\n ", &numel, &numnp, &nmat, &nmode, &plane_stress_flag);
        fprintf( o2, "    %4d %4d %4d %4d %4d\n ", numel, numnp, nmat, nmode, plane_stress_flag);
        fprintf( o3, "    %4d %4d %4d %4d %4d\n ", numel, numnp, nmat, nmode, plane_stress_flag);
        fprintf( o4, "    %4d %4d %4d %4d %4d\n ", numel, numnp, nmat, nmode, plane_stress_flag);
        fgets( buf, BUFSIZ, o1 );

        fprintf( o2, "matl no., E modulus, Poisson Ratio, density \n");
        fprintf( o3, "matl no., E modulus, Poisson Ratio, density \n");
        fprintf( o4, "matl no., E modulus, Poisson Ratio, density \n");
        for( i = 0; i < nmat; ++i )
        {
           fscanf( o1, "%d\n ",&dum);
           fprintf( o2, "%3d ",dum);
           fprintf( o3, "%3d ",dum);
           fprintf( o4, "%3d ",dum);
           fscanf( o1, " %lf %lf %lf\n",&fdum1, &fdum2, &fdum3);
           fprintf( o2, " %9.4f %9.4f %9.4f\n ",fdum1, fdum2, fdum3);
           fprintf( o3, " %9.4f %9.4f %9.4f\n ",fdum1, fdum2, fdum3);
           fprintf( o4, " %9.4f %9.4f %9.4f\n ",fdum1, fdum2, fdum3);
        }
        fgets( buf, BUFSIZ, o1 );

        fprintf( o2, "el no., connectivity, matl no. \n");
        fprintf( o3, "el no., connectivity, matl no. \n");
        fprintf( o4, "el no., connectivity, matl no. \n");
        for( i = 0; i < numel; ++i )
        {
           fscanf( o1,"%d ",&dum);
           fprintf( o2, "%4d ",dum);
           fprintf( o3, "%4d ",dum);
           fprintf( o4, "%4d ",dum);
           for( j = 0; j < npel; ++j )
           {
                fscanf( o1, "%d",&dum3);
                fprintf( o2, "%4d ",dum3);
                fprintf( o3, "%4d ",dum3);
                fprintf( o4, "%4d ",dum3);
           }
           fscanf( o1,"%d\n",&dum4);
           fprintf( o2, " %3d\n",dum4);
           fprintf( o3, " %3d\n",dum4);
           fprintf( o4, " %3d\n",dum4);
        }
        fgets( buf, BUFSIZ, o1 );

        fprintf( o2, "node no., coordinates \n");
        fprintf( o3, "node no., coordinates \n");
        fprintf( o4, "node no., coordinates \n");
        for( i = 0; i < numnp; ++i )
        {
           fscanf( o1,"%d ",&dum);
           fprintf( o2, "%d ",dum);
           fprintf( o3, "%d ",dum);
           fprintf( o4, "%d ",dum);
           fscanf( o1, "%lf %lf",&fdum1,&fdum2);
	   fdum3 = .6923*fdum1 - .80333*fdum2;
           fprintf( o2, "%9.4f %9.4f %9.4f",0.0,fdum2,fdum3);
           fprintf( o3, "%9.4f %9.4f %9.4f",fdum1,0.0,fdum3);
           fprintf( o4, "%9.4f %9.4f %9.4f",fdum1,fdum2,0.0);

/* For patch.x */
           *(ux2+i) =  0.0;
           *(uy2+i) = -0.020*fdum2 + 0.035*fdum3;
           *(uz2+i) =  0.033*fdum2 - 0.027*fdum3;

/* For patch.y */
	   *(ux3+i) = 0.003*fdum1 + 0.021*fdum3;
	   *(uy3+i) = 0.0;
	   *(uz3+i) = 0.041*fdum1 - 0.027*fdum3;

/* For patch.z */
	   *(ux4+i) = 1.000e-04*fdum1 + 3.000e-04*fdum2;
	   *(uy4+i) = 2.000e-04*fdum1 + 4.000e-04*fdum2;
	   *(uz4+i) = 0.0;

           fscanf( o1,"\n");
           fprintf( o2, "\n");
           fprintf( o3, "\n");
           fprintf( o4, "\n");
        }
        fgets( buf, BUFSIZ, o1 );

        dum= 0;
        fprintf( o2, "prescribed displacement x: node  disp value\n");
        fprintf( o3, "prescribed displacement x: node  disp value\n");
        fprintf( o4, "prescribed displacement x: node  disp value\n");
        for( i = 0; i < 4; ++i )
        {
                fprintf( o2, "%4d %14.6e\n",i,*(ux2+i));
                fprintf( o3, "%4d %14.6e\n",i,*(ux3+i));
                fprintf( o4, "%4d %14.6e\n",i,*(ux4+i));
        }
        fprintf( o2, "%4d %14.6e\n",i,0.0);
        for( i = 5; i < numnp; ++i )
        {
                fprintf( o2, "%4d %14.6e\n",i,*(ux2+i));
                fprintf( o3, "%4d %14.6e\n",i,*(ux3+i));
                fprintf( o4, "%4d %14.6e\n",i,*(ux4+i));
        }
        fprintf( o2, "%4d\n ",-10);
        fprintf( o3, "%4d\n ",-10);
        fprintf( o4, "%4d\n ",-10);
        fprintf( o2, "prescribed displacement y: node  disp value\n");
        fprintf( o3, "prescribed displacement y: node  disp value\n");
        fprintf( o4, "prescribed displacement y: node  disp value\n");
        for( i = 0; i < 4; ++i )
        {
                fprintf( o2, "%4d %14.6e\n",i,*(uy2+i));
                fprintf( o3, "%4d %14.6e\n",i,*(uy3+i));
                fprintf( o4, "%4d %14.6e\n",i,*(uy4+i));
        }
        fprintf( o3, "%4d %14.6e\n",i,0.0);
        for( i = 5; i < numnp; ++i )
        {
                fprintf( o2, "%4d %14.6e\n",i,*(uy2+i));
                fprintf( o3, "%4d %14.6e\n",i,*(uy3+i));
                fprintf( o4, "%4d %14.6e\n",i,*(uy4+i));
        }
        fprintf( o2, "%4d\n ",-10);
        fprintf( o3, "%4d\n ",-10);
        fprintf( o4, "%4d\n ",-10);
        fprintf( o2, "prescribed displacement z: node  disp value\n");
        fprintf( o3, "prescribed displacement z: node  disp value\n");
        fprintf( o4, "prescribed displacement z: node  disp value\n");
        for( i = 0; i < 4; ++i )
        {
                fprintf( o2, "%4d %14.6e\n",i,*(uz2+i));
                fprintf( o3, "%4d %14.6e\n",i,*(uz3+i));
                fprintf( o4, "%4d %14.6e\n",i,*(uz4+i));
        }
        fprintf( o4, "%4d %14.6e\n",i,0.0);
        for( i = 5; i < numnp; ++i )
        {
                fprintf( o2, "%4d %14.6e\n",i,*(uz2+i));
                fprintf( o3, "%4d %14.6e\n",i,*(uz3+i));
                fprintf( o4, "%4d %14.6e\n",i,*(uz4+i));
        }
        fprintf( o2, "%4d\n ",-10);
        fprintf( o3, "%4d\n ",-10);
        fprintf( o4, "%4d\n ",-10);

        fgets( buf, BUFSIZ, o1 );
        fgets( buf, BUFSIZ, o1 );
        fprintf( o2, "node with point load and load vector in x,y,z\n");
        fprintf( o3, "node with point load and load vector in x,y,z\n");
        fprintf( o4, "node with point load and load vector in x,y,z\n");
        fprintf( o2, "%4d\n ",-10);
        fprintf( o3, "%4d\n ",-10);
        fprintf( o4, "%4d\n ",-10);
        fprintf( o2, "element and gauss pt. with stress and stress vector in xx,yy,xy\n");
        fprintf( o3, "element and gauss pt. with stress and stress vector in xx,yy,xy\n");
        fprintf( o4, "element and gauss pt. with stress and stress vector in xx,yy,xy\n");
        fprintf( o2, "%4d ",-10);
        fprintf( o3, "%4d ",-10);
        fprintf( o4, "%4d ",-10);

        fprintf( o2, "\n\n\n Note: Remove the prescribed displacements for node 4.\n");
        fprintf( o3, "\n\n\n Note: Remove the prescribed displacements for node 4.\n");
        fprintf( o4, "\n\n\n Note: Remove the prescribed displacements for node 4.\n");
        fprintf( o2, "%4d %14.6e\n",4,*(ux2+4));
        fprintf( o2, "%4d %14.6e\n",4,*(uy2+4));
        fprintf( o2, "%4d %14.6e\n",4,*(uz2+4));
        fprintf( o3, "%4d %14.6e\n",4,*(ux3+4));
        fprintf( o3, "%4d %14.6e\n",4,*(uy3+4));
        fprintf( o3, "%4d %14.6e\n",4,*(uz3+4));
        fprintf( o4, "%4d %14.6e\n",4,*(ux4+4));
        fprintf( o4, "%4d %14.6e\n",4,*(uy4+4));
        fprintf( o4, "%4d %14.6e\n",4,*(uz4+4));

        return 1;
}

