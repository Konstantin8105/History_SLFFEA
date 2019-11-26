/*
    This program performs finite element analysis
    on a brick by reading in the data, assembling the P vector, 
    and then solving the system using dynamic relaxation.  This
    element is capable of handling large deformation loads.
    The material is assumed to be hypo-elastic and the stress
    is updated using the Jaumann Stress Rate.

		Udated 6/20/02

    SLFFEA source file
    Version:  1.2
    Copyright (C) 1999, 2000, 2001  San Le 

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../brick/brconst.h"
#include "../brick/brstruct.h"

#define flag        3                     /* flag for postprocessor */

int brConnectSurfwriter ( int *, int *, char *);

int brVolume( int *, double *, double *);

int brwriter ( BOUND , int *, double *, int *, double *, int *, MATL *,
        char *, STRAIN *, SDIM *, STRESS *, SDIM *, double *);

int matX(double *,double *,double *, int ,int ,int );

int brPassemble(int *, double *, double *, int *, MATL *, double *, STRESS *,
	double * );

int brFMassemble(int *, double *, double *, int *, double *, double *, MATL *,
	double *);

int brformid( BOUND, int *);

int brshg( double *, int, double *, double *, double *);

int brreader( BOUND , int *, double *, int *, double *, MATL *, char *, FILE *,
        STRESS *, SDIM *, double *);

int brMemory( double **, int, int **, int, MATL **, int , XYZI **, int,
        SDIM **, int, STRAIN **, STRESS **, int );

int brshl(double, double *, double * );

int dof, neqn, nmat, nmode, numnp, numel, numel_surf, sdof, sof;
int stress_read_flag, element_stress_read_flag, element_stress_print_flag,
	gauss_stress_flag;

int standard_flag;
       
double shg[sosh], shgh[sosh], shl[sosh],w[num_int], *Vol0;
double dt,cc;
int iteration;
double tolerance;

int main(int argc, char** argv)
{
	int i, j, k, i2;
	int *id, check, iteration_max, name_length, counter, MemoryCounter, dum;
	double  d_iteration_max, iteration_const;
	int matl_num;
        double Emod, G, K, Pois;
	XYZI *mem_XYZI;
	int *mem_int, sofmi, sofmf, sofmSTRESS, sofmXYZI, sofmSDIM, ptr_inc;
	MATL *matl;
	double *mem_double;
        double fpointx, fpointy, fpointz;
	int *connect, *connect_surf, *el_matl, *el_matl_surf, node;
	double *coord, *coord0, *force,*P_global,*U, *dU, *V, *dUm1,
		*Voln, *mass, *node_counter;
        double coord_el_trans[neqel], *coordh;
	char name[30], buf[ BUFSIZ ];
        FILE *o1, *o2, *o3, *o4;
	double det[num_int], wXdet;
	BOUND bc;
	STRESS *stress;
	STRAIN *strain;
	SDIM *stress_node, *strain_node, *mem_SDIM;
        double invariant[nsd], yzsq, zxsq, xysq, xxyy;
	double g;
	double hydro;
	long timec;
	long timef;
	int loc;
	double v,bet,alp,timer;
	double RAM_max, RAM_usable, RAM_needed, MEGS;
	int surf_el_counter, surface_el_flag;

        sof = sizeof(double);

/* Create local shape funcions */

	g = 2.0/sq3;
        check = brshl(g, shl, w );
        if(!check) printf( "Problems with brshl \n");

        memset(name,0,30*sizeof(char));
	
    	printf("What is the name of the file containing the \n");
    	printf("brick structural data? \n");
    	scanf( "%30s",name);

/*   o1 contains all the structural data */
/*   o2 contains input parameters */
/*   o4 contains dynamic data   */

        o1 = fopen( name,"r" );
        o2 = fopen( "brinput","r" );
        o3 = fopen( "data","w" );
        o4 = fopen( "brdynam.dat","r" );

	if(o1 == NULL ) {
		printf("Can't find file %30s\n",name);
		exit(1);
	}

	if( o2 == NULL ) {
		printf("Can't find file brinput\n");
		tolerance = 1.e-13;
		iteration_max = 2000;
		RAM_max = 160.0;
		element_stress_read_flag = 0;
		element_stress_print_flag = 0;
		gauss_stress_flag = 1;
	}
	else
	{
		fgets( buf, BUFSIZ, o2 );
		fscanf( o2, "%lf\n ",&tolerance);
		fgets( buf, BUFSIZ, o2 );
		fscanf( o2, "%d\n ",&iteration_max);
		fgets( buf, BUFSIZ, o2 );
		fscanf( o2, "%lf\n ",&RAM_max);
		fgets( buf, BUFSIZ, o2 );
		fscanf( o2, "%d\n ",&element_stress_read_flag);
		fgets( buf, BUFSIZ, o2 );
		fscanf( o2, "%d\n ",&element_stress_print_flag);
		fgets( buf, BUFSIZ, o2 );
		fscanf( o2, "%d\n ",&gauss_stress_flag);
	}
	gauss_stress_flag = 1;

	if( o4 == NULL ) {
		printf("Can't find file brdynam.dat\n");
		exit(1);
	}

        fgets( buf, BUFSIZ, o1 );
        fscanf( o1, "%d %d %d %d\n ",&numel,&numnp,&nmat,&nmode);
        dof=numnp*ndof;
        sdof=numnp*nsd;
	standard_flag = 1;

/*   Begin allocation of meomory */

	MemoryCounter = 0;

/* For the doubles */
	sofmf=3*numnp*nsd + 7*dof + 2*numel + numnp;
	MemoryCounter += sofmf*sizeof(double);
	printf( "\n Memory requrement for doubles is %15d bytes\n",MemoryCounter);

/* For the integers */
	sofmi=2*numel*npel + dof + 2*numel + numnp+1;
	MemoryCounter += sofmi*sizeof(int);
	printf( "\n Memory requrement for integers is %15d bytes\n",MemoryCounter);

/* For the XYZI integers */
	sofmXYZI=numnp+1;
	MemoryCounter += sofmXYZI*sizeof(XYZI);
	printf( "\n Memory requrement for XYZI integers is %15d bytes\n",MemoryCounter);

/* For the SDIM doubles */
	sofmSDIM = 2*numnp;
	MemoryCounter += sofmSDIM*sizeof(SDIM);
	printf( "\n Memory requrement for SDIM doubles is %15d bytes\n",MemoryCounter);

/* For the STRESS */
	sofmSTRESS=numel;
	MemoryCounter += sofmSTRESS*sizeof(STRESS) + sofmSTRESS*sizeof(STRAIN);
	printf( "\n Memory requrement for STRESS doubles is %15d bytes\n",MemoryCounter);

	check = brMemory( &mem_double, sofmf, &mem_int, sofmi, &matl, nmat,
		&mem_XYZI, sofmXYZI, &mem_SDIM, sofmSDIM, &strain, &stress,
		sofmSTRESS );

	if(!check) printf( "Problems with brMemory \n");

/* For the doubles */
					        ptr_inc=0; 
	coord0=(mem_double+ptr_inc);            ptr_inc += numnp*nsd;
	coord=(mem_double+ptr_inc);             ptr_inc += numnp*nsd;
	coordh=(mem_double+ptr_inc);            ptr_inc += numnp*nsd;
	force=(mem_double+ptr_inc);             ptr_inc += dof;
	U=(mem_double+ptr_inc);                 ptr_inc += dof;
	dU=(mem_double+ptr_inc);                ptr_inc += dof;
	dUm1=(mem_double+ptr_inc);              ptr_inc += dof;
	V=(mem_double+ptr_inc);                 ptr_inc += dof;
	mass=(mem_double+ptr_inc);              ptr_inc += dof;
	P_global=(mem_double+ptr_inc); 	        ptr_inc += dof;
	Voln=(mem_double+ptr_inc);             	ptr_inc += numel; 
	Vol0=(mem_double+ptr_inc);              ptr_inc += numel;
	node_counter=(mem_double+ptr_inc);      ptr_inc += numnp; 

/* For the integers */
						ptr_inc = 0; 
	connect=(mem_int+ptr_inc); 		ptr_inc += numel*npel; 
	connect_surf=(mem_int+ptr_inc);         ptr_inc += numel*npel;
	id=(mem_int+ptr_inc);                   ptr_inc += dof;
        el_matl=(mem_int+ptr_inc);              ptr_inc += numel;
	el_matl_surf=(mem_int+ptr_inc);         ptr_inc += numel;
        bc.force =(mem_int+ptr_inc);            ptr_inc += numnp;
        bc.num_force=(mem_int+ptr_inc);         ptr_inc += 1;

/* For the XYZI integers */
					  ptr_inc = 0; 
	bc.fix =(mem_XYZI+ptr_inc);       ptr_inc += numnp;
	bc.num_fix=(mem_XYZI+ptr_inc);    ptr_inc += 1;

/* For the SDIM doubles */
                                                ptr_inc = 0;
	stress_node=(mem_SDIM+ptr_inc);         ptr_inc += numnp;
	strain_node=(mem_SDIM+ptr_inc);         ptr_inc += numnp;

	timec = clock();
	timef = 0;

        for( k = 0; k < numel; ++k )
        {
        	for( j = 0; j < num_int; ++j )
        	{
                 	stress[k].pt[j].xx=.0;
                 	stress[k].pt[j].yy=.0;
                 	stress[k].pt[j].zz=.0;
                 	stress[k].pt[j].xy=.0;
                 	stress[k].pt[j].zx=.0;
                 	stress[k].pt[j].yz=.0;
		}
	}

	stress_read_flag = 1;
        check = brreader( bc, connect, coord0, el_matl, force, matl, name, o1,
        	stress, stress_node, U);
	if(!check) printf( "Problems with reader \n");

	printf(" \n\n");

	printf( "\n The total memory requrement is %15d bytes\n", MemoryCounter);
	MEGS = ((double)(MemoryCounter))/MB;
	printf( "\n Which is %16.4e MB\n", MEGS);

/* Because this is a dynamic code, the id array is not needed except for
   the brwriter subroutine */

	check = brformid( bc, id );
	if(!check) printf( "Problems with brformid \n");

	for( i = 0; i < dof; ++i )
	{
/* Initializing the data for change in displacement */

		 *(dU+i)=.0;
		 *(dUm1+i)=.0;

/* Initializing the data for velocity */

		*(V+i)=0.0;
	}

	for( i = 0; i < sdof; ++i )
	{

/* Initializing the data for x,y,z coordindates */

		*(coord+i)=*(coord0+i);
		*(coordh+i)=*(coord0+i);
	}

/* Initializing the data for the Volume */

	check = brVolume( connect, coord0, Vol0);
        if(!check) printf( "Problems with brVolume \n");

	for( k = 0; k < numel; ++k )
	{
       	   	printf(" \n");
	   	printf("Volume for element %d %14.5f ", k, *(Vol0+k));
	}
	printf("\n\n");

/* The mass matrix */

       	printf(" \n");
	printf("mass  matrix\n");

    	check = brFMassemble(connect, coord, coordh, el_matl, force, mass, matl, U );
	if(!check) printf( "Problems with brFMassemble \n");

	printf( "\n\n This is the mass matrix \n");
	for( i = 0; i < dof; ++i )
	{
	    printf( "\ndof %5d   %14.5f",i,*(mass+i));
	}
	printf(" \n");

	printf(" \n");

	fscanf( o4, " %lf \n", &dt);
	fscanf( o4, " %lf \n", &cc);
	fscanf( o4, " %d \n", &iteration_max);

	v=cc*dt;
	bet=(2.0 - v)/(2.0 + v);
	alp=2.0*dt*dt/(2.0 + v);
	timer = 0.;
	d_iteration_max = (double) iteration_max;
	iteration_const = 1.0/(d_iteration_max*dt);

/* Dynamic Relaxation occurs below */

	fprintf(o3, " %d %d \n",iteration_max,1);
	for( iteration = 0; iteration < iteration_max; ++iteration )
	{
		timer=timer+dt;
		/*fprintf(o3, " %f %f \n",timer,*(dU+dof-1)/numel*10);*/
		printf( " \n %3d \n",iteration);
		check = brPassemble( connect, coord, coordh, el_matl, matl, P_global,
		    stress, V);
		if(!check) printf( "Problems with brPassemble \n");
		/*for( j = 0; j < dof; ++j )
		{
		    printf("\n vvvv  %3d %12.10f  %12.10f %14.5f",j,*(dU+j),
			*(V+j),*(P_global+j));
		}*/

		for( j = 0; j < numnp; ++j )
		{
		    for( i2 = 0; i2 < ndof; ++i2 )
		    {
			/*printf("\n  %3d %12.10f  %12.10f %14.5f",
				ndof*j+i2,*(dU+ndof*j+i2),
				*(V+ndof*j),*(P_global+ndof*j+i2));*/
			/*printf( "%4d %f %f %f %f\n",ndof*j+i2,
				*(P_global+ndof*j+i2),*(dU+ndof*j+i2),
				*(V+ndof*j),*(coord+nsd*j));*/
			loc = ndof*j+i2;
			*(dU+loc) = 0.0;
			*(dU+loc) += bet*(*(dU+loc) - *(dUm1+loc)) -
				alp*( *(P_global+loc) - *(force+loc) )/(*(mass+loc));
		    }

/* Update of the XYZ coordinate matrix  */

		    *(coordh+nsd*j) = *(coord+nsd*j)+.5*(*(dU+ndof*j));
		    *(coordh+nsd*j+1) = *(coord+nsd*j+1)+.5*(*(dU+ndof*j+1));
		    *(coordh+nsd*j+2) = *(coord+nsd*j+2)+.5*(*(dU+ndof*j+2));

		    *(coord+nsd*j) += *(dU+ndof*j);
		    *(coord+nsd*j+1) += *(dU+ndof*j+1);
		    *(coord+nsd*j+2) += *(dU+ndof*j+2);

/* Update of the Velocity matrix  */

		    *(V+ndof*j) = *(dU+ndof*j)/dt;
		    *(V+ndof*j+1) = *(dU+ndof*j+1)/dt;
		    *(V+ndof*j+2) = *(dU+ndof*j+2)/dt;

/* Set dUm1 to -dU  */

		    *(dUm1+ndof*j)=-*(dU+ndof*j);
		    *(dUm1+ndof*j+1)=-*(dU+ndof*j+1);
		    *(dUm1+ndof*j+2)=-*(dU+ndof*j+2);

		}
/*
   To handle prescribed displacements, I add each displacement
   incrementally to their corresponding nodes based on the following:

    timer = iteration*dt

         V = U_fixed/(iteration_max*dt)
     coord = coord0 + [U_fixed/(iteration_max*dt)]*timer
    coordh = coord0 + [U_fixed/(iteration_max*dt)]*(timer - dt/2.0)
*/
		for( i2 = 0; i2 < bc.num_fix[0].x; ++i2 )
		{
		    node=bc.fix[i2].x;
		    *(V+ndof*node)=*(U+ndof*node)*iteration_const;
		    *(coord+nsd*node)=*(coord0+nsd*node) +
			*(U+ndof*node)*iteration_const*timer;
		    *(coordh+nsd*node)=*(coord0+nsd*node) +
			*(U+ndof*node)*iteration_const*(timer - dt/2.0);
		}
		for( i2 = 0; i2 < bc.num_fix[0].y; ++i2 )
		{
		    node=bc.fix[i2].y;
		    *(V+ndof*node+1)=*(U+ndof*node+1)*iteration_const;
		    *(coord+nsd*node+1)=*(coord0+nsd*node+1) +
			*(U+ndof*node+1)*iteration_const*timer;
		    *(coordh+nsd*node+1)=*(coord0+nsd*node+1) +
			*(U+ndof*node+1)*iteration_const*(timer - dt/2.0);
		}
		for( i2 = 0; i2 < bc.num_fix[0].z; ++i2 )
		{
		    node=bc.fix[i2].z;
		    *(V+ndof*node+2)=*(U+ndof*node+2)*iteration_const;
		    *(coord+nsd*node+2)=*(coord0+nsd*node+2) +
			*(U+ndof*node+2)*iteration_const*timer;
		    *(coordh+nsd*node+2)=*(coord0+nsd*node+2) +
			*(U+ndof*node+2)*iteration_const*(timer - dt/2.0);
		    /*printf("\n vvvv  %3d %3d %12.10f  %12.10f %14.5f %14.5f",i2,iteration,
			*(coord+nsd*node+2), *(coordh+nsd*node+2), *(U+ndof*node+2),
			iteration_const);*/
		}
	}


/* Calculating the final value of the Volumes */

	check = brVolume( connect, coord, Voln);
        if(!check) printf( "Problems with brVolume \n");

/*
        printf("\nThis is the Volume\n");
        for( i = 0; i < numel; ++i )
        {
                printf("%4i %12.4e\n",i, *(Voln + i));
        }
*/

/* Set reaction forces */

	for( i = 0; i < numnp; ++i )
	{
	   if( *(id+ndof*i) < 0 || *(id+ndof*i+1) < 0 || *(id+ndof*i+2) < 0 )
           {
		*(force+ndof*i) = *(P_global+ndof*i);
		*(force+ndof*i+1) = *(P_global+ndof*i+1);
		*(force+ndof*i+2) = *(P_global+ndof*i+2);
	   }
	}

        printf(" \n");
        /*printf( "%f %f %f \n",*(force+7),*(P_global+7),*(dU+7));*/
       	printf(" \n\n");
        for( i = 0; i < numnp; i+=4 )
        {
		printf("\n %3d        %7.4f    %7.4f    %7.4f",i+1,
			*(coord+nsd*i), *(coord+nsd*i+1),*(coord+nsd*i+2));
		printf("\n %3d        %7.4f    %7.4f    %7.4f",i+2,
			*(coord+nsd*(i+3)), *(coord+nsd*(i+3)+1),*(coord+nsd*(i+3)+2));
		printf("\n %3d        %7.4f    %7.4f    %7.4f", i+3,
			*(coord+nsd*(i+1)), *(coord+nsd*(i+1)+1),*(coord+nsd*(i+1)+2));
		printf("\n %3d        %7.4f    %7.4f    %7.4f",i+4,
			*(coord+nsd*(i+2)), *(coord+nsd*(i+2)+1),*(coord+nsd*(i+2)+2));
	}
       	printf(" \n\n");
        for( k = 0; k < numel; ++k )
        {
                matl_num = *(el_matl+k);
                Emod = matl[matl_num].E;
                Pois = matl[matl_num].nu;

                K=Emod/(1.0-2*Pois)/3.0;
                G=Emod/(1.0+Pois)/2.0;

        	for( j = 0; j < num_int; ++j )
        	{

/* Count the number of times a particular node is part of an element */

		    node = *(connect+npel*k+j);

                    *(node_counter + node) += 1.0;

		    hydro = (stress[k].pt[j].xx + stress[k].pt[j].yy +
			stress[k].pt[j].zz)*Pois/Emod;
		    strain[k].pt[j].xx = stress[k].pt[j].xx/(2.0*G) - hydro;
                    strain[k].pt[j].yy = stress[k].pt[j].yy/(2.0*G) - hydro;
                    strain[k].pt[j].zz = stress[k].pt[j].zz/(2.0*G) - hydro;
		    strain[k].pt[j].xy = stress[k].pt[j].xy/G;
                    strain[k].pt[j].zx = stress[k].pt[j].zx/G;
                    strain[k].pt[j].yz = stress[k].pt[j].yz/G;

/* Calculate the principal strians */

                    memset(invariant,0,nsd*sof);
                    xysq = strain[k].pt[j].xy*strain[k].pt[j].xy;
                    zxsq = strain[k].pt[j].zx*strain[k].pt[j].zx;
                    yzsq = strain[k].pt[j].yz*strain[k].pt[j].yz;
                    xxyy = strain[k].pt[j].xx*strain[k].pt[j].yy;

                    *(invariant) = - strain[k].pt[j].xx -
                         strain[k].pt[j].yy - strain[k].pt[j].zz;
                    *(invariant+1) = xxyy +
                         strain[k].pt[j].yy*strain[k].pt[j].zz +
                         strain[k].pt[j].zz*strain[k].pt[j].xx - yzsq - zxsq - xysq;
                    *(invariant+2) =
                         - xxyy*strain[k].pt[j].zz -
                         2*strain[k].pt[j].yz*strain[k].pt[j].zx*strain[k].pt[j].xy +
                         yzsq*strain[k].pt[j].xx + zxsq*strain[k].pt[j].yy +
                         xysq*strain[k].pt[j].zz;

                    check = cubic(invariant);

                    strain[k].pt[j].I = *(invariant);
                    strain[k].pt[j].II = *(invariant+1);
                    strain[k].pt[j].III = *(invariant+2);


/* Add all the strains for a particular node from all the elements which share that node */

		    strain_node[node].xx += strain[k].pt[j].xx;
		    strain_node[node].yy += strain[k].pt[j].yy;
		    strain_node[node].zz += strain[k].pt[j].zz;
		    strain_node[node].xy += strain[k].pt[j].xy;
		    strain_node[node].zx += strain[k].pt[j].zx;
		    strain_node[node].yz += strain[k].pt[j].yz;
		    strain_node[node].I += strain[k].pt[j].I;
		    strain_node[node].II += strain[k].pt[j].II;
		    strain_node[node].III += strain[k].pt[j].III;

/* Calculate the principal stresses */

                    memset(invariant,0,nsd*sof);
                    xysq = stress[k].pt[j].xy*stress[k].pt[j].xy;
                    zxsq = stress[k].pt[j].zx*stress[k].pt[j].zx;
                    yzsq = stress[k].pt[j].yz*stress[k].pt[j].yz;
                    xxyy = stress[k].pt[j].xx*stress[k].pt[j].yy;

                    *(invariant) = - stress[k].pt[j].xx -
                         stress[k].pt[j].yy - stress[k].pt[j].zz;
                    *(invariant+1) = xxyy +
                         stress[k].pt[j].yy*stress[k].pt[j].zz +
                         stress[k].pt[j].zz*stress[k].pt[j].xx - yzsq - zxsq - xysq;
                    *(invariant+2) =
                         - xxyy*stress[k].pt[j].zz -
                         2*stress[k].pt[j].yz*stress[k].pt[j].zx*stress[k].pt[j].xy +
                         yzsq*stress[k].pt[j].xx + zxsq*stress[k].pt[j].yy +
                         xysq*stress[k].pt[j].zz;

                    check = cubic(invariant);

                    stress[k].pt[j].I = *(invariant);
                    stress[k].pt[j].II = *(invariant+1);
                    stress[k].pt[j].III = *(invariant+2);

/* Add all the stresses for a particular node from all the elements which share that node */

		    stress_node[node].xx += stress[k].pt[j].xx;
		    stress_node[node].yy += stress[k].pt[j].yy;
		    stress_node[node].zz += stress[k].pt[j].zz;
		    stress_node[node].xy += stress[k].pt[j].xy;
		    stress_node[node].zx += stress[k].pt[j].zx;
		    stress_node[node].yz += stress[k].pt[j].yz;
		    stress_node[node].I += stress[k].pt[j].I;
		    stress_node[node].II += stress[k].pt[j].II;
		    stress_node[node].III += stress[k].pt[j].III;
		}
	}

/* Average all the stresses and strains at the nodes */

	for( i = 0; i < numnp ; ++i )
	{
		strain_node[i].xx /= *(node_counter + i);
		strain_node[i].yy /= *(node_counter + i);
		strain_node[i].zz /= *(node_counter + i);
		strain_node[i].xy /= *(node_counter + i);
		strain_node[i].zx /= *(node_counter + i);
		strain_node[i].yz /= *(node_counter + i);
		strain_node[i].I /= *(node_counter + i);
		strain_node[i].II /= *(node_counter + i);
		strain_node[i].III /= *(node_counter + i);

		stress_node[i].xx /= *(node_counter + i);
		stress_node[i].yy /= *(node_counter + i);
		stress_node[i].zz /= *(node_counter + i);
		stress_node[i].xy /= *(node_counter + i);
		stress_node[i].zx /= *(node_counter + i);
		stress_node[i].yz /= *(node_counter + i);
		stress_node[i].I /= *(node_counter + i);
		stress_node[i].II /= *(node_counter + i);
		stress_node[i].III /= *(node_counter + i);
	}

        for( k = 0; k < numel; ++k )
        {
		surface_el_flag = 0;
		for( j = 0; j < num_int; ++j )
		{

/* Determine which elements are on the surface */

		   node = *(connect+npel*k+j);

		   if( *(node_counter + node) < 7.5 )
			surface_el_flag = 1;
		}

		if(surface_el_flag)
		{
		    for( j = 0; j < npel; ++j)
		    {
			*(connect_surf + npel*surf_el_counter + j) =
				*(connect + npel*k + j);
		    }
		    *(el_matl_surf + surf_el_counter) = matl_num;
		    ++surf_el_counter;
		}
		numel_surf = surf_el_counter;
	}

        for( i = 0; i < numnp; ++i )
        {
	   if( *(id+ndof*i) > -1 )
		*(U+ndof*i) = *(coord+nsd*i)-*(coord0+nsd*i); 
	   if( *(id+ndof*i+1) > -1 )
		*(U+ndof*i+1) = *(coord+nsd*i+1)-*(coord0+nsd*i+1);
	   if( *(id+ndof*i+2) > -1 )
		*(U+ndof*i+2) = *(coord+nsd*i+2)-*(coord0+nsd*i+2);
	}

        check = brwriter ( bc, connect, coord0, el_matl, force, id, matl,
                name, strain, strain_node, stress, stress_node, U);
        if(!check) printf( "Problems with brwriter \n");

	check = brConnectSurfwriter( connect_surf, el_matl_surf, name);
	if(!check) printf( "Problems with brConnectSurfwriter \n");

        for( i = 0; i < numnp; ++i )
        {
           	printf("\n  %3d   %14.5f  %14.5f  %14.5f",nsd*i,*(coord+nsd*i),
			*(P_global+nsd*i), *(force+nsd*i));
           	printf("\n  %3d   %14.5f  %14.5f  %14.5f",nsd*i+1,*(coord+nsd*i+1),
			*(P_global+nsd*i+1), *(force+nsd*i+1));
        	printf("\n  %3d   %14.5f  %14.5f  %14.5f",nsd*i+2,*(coord+nsd*i+2),
			*(P_global+nsd*i+2), *(force+nsd*i+2));
	}
       	printf("\n");

	timec = clock();
	printf("elapsed CPU = %f\n",( (double)timec)/800.);

	free(strain);
	free(stress);
	free(matl);
	free(mem_double);
	free(mem_int);
	free(mem_XYZI);
}
