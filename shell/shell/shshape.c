/*
     SLFFEA source file
     Version:  1.3
     Copyright (C) 1999, 2000, 2001, 2002  San Le

     The source code contained in this file is released under the
     terms of the GNU Library General Public License.
*/

#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "shconst.h"
#include "shstruct.h"

int normal(double *);

int normcrossX(double *, double *, double *);

int dotX(double *,double *, double *, int );

int shshl( double g, SH shl, double *w)
{
/* 
     This subroutine calculates the local shape function derivatives for
     a shell element at the gauss points.

     It is based on the subroutine QDCSHL from the book
     "The Finite Element Method" by Thomas Hughes, page 784.

			Updated 5/5/99
*/
        double ra[]={-0.50, 0.50, 0.50,-0.50, 0.00};
        double sa[]={-0.50,-0.50, 0.50, 0.50, 0.00};
        double ta[]={-0.50, 0.50 };
	double temps,tempr,tempt,r,s,t;
	int i,j,k;

        for( k = 0; k < num_intb; ++k )
        {
/* 
   calculating the weights and local dN/ds,dr matrix for 2X2 
   integration for shl
*/
		*(w+k)=1.0;
		*(w+num_intb+k)=1.0;

        	r=g*(*(ra+k));
        	s=g*(*(sa+k));
        	for( i = 0; i < npell; ++i )
        	{
		   tempr = pt5 + *(ra+i)*r;
		   temps = pt5 + *(sa+i)*s;
        	   shl.bend[npell*(nsdl+1)*k+i]=*(ra+i)*temps;
        	   shl.bend[npell*(nsdl+1)*k+npell*1+i]=tempr*(*(sa+i));
        	   shl.bend[npell*(nsdl+1)*k+npell*2+i]=tempr*temps;
		}
               	/*printf("\n");*/
	}
/* 
   calculating the weights and local dN/ds,dr matrix for 1X1 
   point integration for shl in lamina
*/
        r=g*(*(ra+num_intb));
        s=g*(*(sa+num_intb));
	/* Actually, 4.0 but I'm adding 4 times in shKassemble */
        for( i = 0; i < npell; ++i )
        {
	   tempr = pt5 + *(ra+i)*r;
	   temps = pt5 + *(sa+i)*s;
           shl.shear[i]=*(ra+i)*temps;
           shl.shear[npell*1+i]=tempr*(*(sa+i));
           shl.shear[npell*2+i]=tempr*temps;
	}

/* 
   calculating the local N and dN/dt matrix for 2X2 point
   integration for shl in fiber 
*/

   t=.25*g;

   shl.bend_z[0]= *(ta); 	         shl.bend_z[1]= *(ta+1);
   shl.bend_z[2]=pt5+t;     	         shl.bend_z[3]=pt5-t;
   shl.bend_z[npelf*num_ints]= *(ta);    shl.bend_z[npelf*num_ints+1]= *(ta+1);
   shl.bend_z[npelf*num_ints+2]=pt5-t;   shl.bend_z[npelf*num_ints+3]=pt5+t;

/* 
   calculating the local N and dN/dt matrix for 1X1 point
   integration for shl in fiber 
*/
   t=.25*g;
   shl.shear_z[0]= *(ta); 	         shl.shear_z[1]= *(ta+1);
   shl.shear_z[npelf*(num_ints-1)]=pt5;  shl.shear_z[npelf*(num_ints-1)+1]=pt5;

        return 1;
}

int shshg( double *det, int el, SH shl, SH shg, XL xl, double *zp1, double *zm1, 
	double *znode, double *dzdt_node, ROTATE rotate)

{
/*
     This subroutine calculates the global shape function derivatives for
     a shell element at the gauss points.

     It is based on the subroutine QDCSHG from the book
     "The Finite Element Method" by Thomas Hughes, page 783.

 ....  CALCULATE GLOBAL DERIVATIVES OF SHAPE FUNCTIONS AND
       JACOBIAN DETERMINANTS FOR A FOUR-NODE QUADRALATERAL ELEMENT

       xl.bar[j+npell*i] = GLOBAL COORDINATES ON LAMINA 
       xl.hat[j+npell*i] = GLOBAL COORDINATES ON FIBER
       *(det+num_intb*k2+k)  = JACOBIAN DETERMINANT
       *(det+num_int+k2)  = JACOBIAN DETERMINANT FOR 1X1 GAUSS

       FOR 2X2 PT. GAUSS
       shl.bend[npell*(nsdl+1)*k+i] = LOCAL ("XI") DERIVATIVE OF SHAPE FUNCTION
       shl.bend[npell*(nsdl+1)*k+npell*1+i] = LOCAL ("ETA") DERIVATIVE OF SHAPE FUNCTION
       shl.bend[npell*(nsdl+1)*k+npell*2+i] = LOCAL SHAPE FUNCTION(in XI and ETA)
       shg.bend[npell*(nsd+1)*k+i] = X-DERIVATIVE OF SHAPE FUNCTION
       shg.bend[npell*(nsd+1)*k+npell*1+i] = Y-DERIVATIVE OF SHAPE FUNCTION
       shg.bend[npell*(nsd+1)*k+npell*2+i] = Z-DERIVATIVE OF SHAPE FUNCTION
       shg.bend[npell*(nsd+1)*k+npell*3+i] = shl(npell*(nsd)*k+npell*2+i)
       shg.bend_z[npell*(nsd)*k+i] = X-DERIVATIVE OF SHAPE FUNCTION
       shg.bend_z[npell*(nsd)*k+npell*1+i] = Y-DERIVATIVE OF SHAPE FUNCTION
       shg.bend_z[npell*(nsd)*k+npell*2+i] = Z-DERIVATIVE OF SHAPE FUNCTION

       FOR 1X1 PT. GAUSS IN LAMINA
       shl.shear[i] = LOCAL ("XI") DERIVATIVE OF SHAPE FUNCTION 
       shl.shear[npell*1+i] = LOCAL ("ETA") DERIVATIVE OF SHAPE FUNCTION
       shl.shear[npell*2+i] = LOCAL SHAPE FUNCTION
       shg.shear[i] = X-DERIVATIVE OF SHAPE FUNCTION
       shg.shear[npell*1+i] = Y-DERIVATIVE OF SHAPE FUNCTION
       shg.shear[npell*2+i] = Z-DERIVATIVE OF SHAPE FUNCTION
       shg.shear[npell*3+i] = shl(npell*(nsd)*num_intb+npell*2+i)
       shg.shear_z[i] = X-DERIVATIVE OF SHAPE FUNCTION
       shg.shear_z[npell*1+i] = Y-DERIVATIVE OF SHAPE FUNCTION
       shg.shear_z[npell*2+i] = Z-DERIVATIVE OF SHAPE FUNCTION

       FOR 2X2 PT. GAUSS
       shl.bend_z[i] = LOCAL ("ZETA") DERIVATIVE OF SHAPE FUNCTION
       shl.bend_z[npelf*num_ints+npelf*k2+i] = LOCAL SHAPE FUNCTION(in ZETA)
       FOR 1X1 AND 2X2 PT. GAUSS  IN LAMINA
       shl.shear_z[i] = LOCAL ("ZETA") DERIVATIVE OF SHAPE FUNCTION
       shl.shear_z[npelf+i] = LOCAL SHAPE FUNCTION(in ZETA) 

       *(xs+2*j+i) = JACOBIAN MATRIX
          i    = LOCAL NODE NUMBER OR GLOBAL COORDINATE NUMBER
          j    = GLOBAL COORDINATE NUMBER
          k    = INTEGRATION-POINT NUMBER FOR LAMINA 
          k2   = INTEGRATION-POINT NUMBER FOR FIBER
       num_intb    = NUMBER OF INTEGRATION POINTS FOR LAMINA, 4
       num_ints    = NUMBER OF INTEGRATION POINTS FOR FIBER, 2
	
			Updated 9/25/01
*/
        double xs[soxshat],temp[nsdsq],col1[nsdl],col2[nsdl],temp1,temp2;
	double shlshl2_vec[npell],xfib2[2];
	double vecl_xi[nsd],vecl_eta[nsd],vecl_alp[nsd],vecl_beta[nsd],
		vecl_1[nsd],vecl_2[nsd],vecl_3[nsd];
	XLXS xlxs;
	int check,i,j,k,k2;

/* 
   calculating the dN/dx,dy,dz dza/dx,dy,dz matrices for 2X2 
   integration
*/
        for( k2 = 0; k2 < num_ints; ++k2 )
	{
           for( i = 0; i < npell; ++i )
	   {
/* 
   calculating the values of the za and dzdt matrix for 2X2 
   integration
*/
	       *(znode+npell*k2+i)=
		    shl.bend_z[npelf*(nsdf+1)*k2+npelf*1+1]*(*(zp1+i))+
		    shl.bend_z[npelf*(nsdf+1)*k2+npelf*1]*(*(zm1+i));
	       *(dzdt_node+i)=
		    shl.bend_z[npelf*(nsdf+1)*k2+1]*(*(zp1+i))+
		    shl.bend_z[npelf*(nsdf+1)*k2]*(*(zm1+i));
	   }
           for( k = 0; k < num_intb; ++k )
	   {
/* 
   calculating dx,dy,dz/ds,dr matrix for 2X2 
   integration
*/
           	for( j = 0; j < nsdl; ++j ) /* j loops over s,r */
	   	{
        		for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
        		{
	     		    check=dotX((xlxs.bar+nsdl*i+j),(shl.bend+npell*(nsdl+1)*k+npell*j),
				(xl.bar+npell*i),npell);
			}
	   	}
               	for( j = 0; j < nsdl; ++j ) /* j loops over s,r */
	       	{
           	   for( i = 0; i < npell; ++i )
	   	   {
	      	   	*(shlshl2_vec+i)=
			   shl.bend[npell*(nsdl+1)*k+npell*j+i]*(*(znode+npell*k2+i));
		   }
              	   for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
        	   {
	     	       check=dotX((xlxs.hat+nsd*i+j),shlshl2_vec,
				(xl.hat+npell*i), npell);
                       *(xs+nsd*i+j)=xlxs.bar[nsdl*i+j]+xlxs.hat[nsd*i+j];
		   }
	       	}
/* 
   calculating dx,dy,dz/dt matrix for 2X2 
   integration
*/
           	for( i = 0; i < npell; ++i )
	   	{
	      	   *(shlshl2_vec+i)=
			shl.bend[npell*(nsdl+1)*k+npell*2+i]*(*(dzdt_node+i));
		}
               	for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
               	{
	           check= dotX((xlxs.hat+nsd*i+2),shlshl2_vec,
			(xl.hat+npell*i), npell);
		   *(xs+nsd*i+2)=xlxs.hat[nsd*i+2];
	       	}
/* 
   calculating rotation matrix for lamina q[i,j] matrix for 2X2 
   integration
*/
		*(vecl_xi)=*(xs);		*(vecl_eta)=*(xs+1);
		*(vecl_xi+1)=*(xs+nsd*1);	*(vecl_eta+1)=*(xs+nsd*1+1);
		*(vecl_xi+2)=*(xs+nsd*2);	*(vecl_eta+2)=*(xs+nsd*2+1);

        	check = normal(vecl_xi);
		if(!check) printf( " Problems with vecl_xi element %d\n",el);
        	check = normal(vecl_eta);
		if(!check) printf( " Problems with vecl_eta element %d\n",el);
        	check = normcrossX(vecl_xi,vecl_eta,vecl_3);
		if(!check) printf( " Problems with normcrossX \n");

		*(vecl_alp)=*(vecl_xi)+*(vecl_eta);
		*(vecl_alp+1)=*(vecl_xi+1)+*(vecl_eta+1);
		*(vecl_alp+2)=*(vecl_xi+2)+*(vecl_eta+2);
        	check = normal(vecl_alp);
		if(!check) printf( " Problems with vecl_alp element %d\n",el);
        	check = normcrossX(vecl_3,vecl_alp,vecl_beta);
		if(!check) printf( " Problems with normcrossX \n");

		*(vecl_1)=rt2*(*(vecl_alp)-*(vecl_beta));
		*(vecl_1+1)=rt2*(*(vecl_alp+1)-*(vecl_beta+1));
		*(vecl_1+2)=rt2*(*(vecl_alp+2)-*(vecl_beta+2));
		*(vecl_2)=rt2*(*(vecl_alp)+*(vecl_beta));
		*(vecl_2+1)=rt2*(*(vecl_alp+1)+*(vecl_beta+1));
		*(vecl_2+2)=rt2*(*(vecl_alp+2)+*(vecl_beta+2));

               	for( i = 0; i < nsd; ++i ) 
               	{
		    rotate.l_bend[nsdsq*num_intb*k2+nsdsq*k+i]=*(vecl_1+i);
		    rotate.l_bend[nsdsq*num_intb*k2+nsdsq*k+nsd*1+i]=*(vecl_2+i);
		    rotate.l_bend[nsdsq*num_intb*k2+nsdsq*k+nsd*2+i]=*(vecl_3+i);
		}

               	*(temp)=*(xs+4)*(*(xs+8))-*(xs+7)*(*(xs+5));
               	*(temp+3)=*(xs+6)*(*(xs+5))-*(xs+3)*(*(xs+8));
               	*(temp+6)=*(xs+3)*(*(xs+7))-*(xs+6)*(*(xs+4));

               	*(det+num_intb*k2+k)=
			*(xs)*(*(temp))+*(xs+1)*(*(temp+3))+*(xs+2)*(*(temp+6));
               	/*printf(" %d %f\n", k, *(det+num_intb*k2+k));*/
 
               	if(*(det+num_intb*k2+k) <= 0.0 )
               	{
                    printf("the element (%d) is inverted; det:%f; integ pt.:%d\n",
                        el,*(det+num_intb*k2+k),k);
                    return 0;
               	}
 
/* The inverse of the jacobian, dc/dx, is calculated below */
 
               	*(temp+1)=*(xs+7)*(*(xs+2))-*(xs+1)*(*(xs+8));
               	*(temp+4)=*(xs)*(*(xs+8))-*(xs+6)*(*(xs+2));
               	*(temp+7)=*(xs+6)*(*(xs+1))-*(xs)*(*(xs+7));
               	*(temp+2)=*(xs+1)*(*(xs+5))-*(xs+4)*(*(xs+2));
               	*(temp+5)=*(xs+3)*(*(xs+2))-*(xs)*(*(xs+5));
               	*(temp+8)=*(xs)*(*(xs+4))-*(xs+3)*(*(xs+1));
     
               	for( j = 0; j < nsd; ++j )
               	{
                    for( i = 0; i < nsd; ++i )
                    {
                       *(xs+nsd*i+j)=*(temp+nsd*i+j)/(*(det+num_intb*k2+k));
                    }
               	}
               	for( i = 0; i < npell; ++i )
               	{
                   *(col1)=shl.bend[npell*(nsd)*k+i];
                   *(col1+1)=shl.bend[npell*(nsd)*k+npell*1+i];
                   *(col2)=*(xs);
                   *(col2+1)=*(xs+nsd*1);
                   check= dotX((shg.bend+npell*(nsd+1)*num_intb*k2+npell*(nsd+1)*k+i),
			(col2),(col1), nsdl);
                   *(col2)=*(xs+1);
                   *(col2+1)=*(xs+nsd*1+1);
                   check= dotX((shg.bend+npell*(nsd+1)*num_intb*k2+npell*(nsd+1)*k+npell*1+i),
			(col2),(col1), nsdl);
                   *(col2)=*(xs+2);
                   *(col2+1)=*(xs+nsd*1+2);
                   check=dotX((shg.bend+npell*(nsd+1)*num_intb*k2+npell*(nsd+1)*k+npell*2+i),
			(col2),(col1), nsdl);
                   shg.bend[npell*(nsd+1)*num_intb*k2+npell*(nsd+1)*k+npell*3+i]=
		      shl.bend[npell*(nsd)*k+npell*2+i];

       		   shg.bend_z[npell*(nsd)*num_intb*k2+npell*(nsd)*k+i] =  
		      *(dzdt_node+i)*(*(xs+nsd*2));
       		   shg.bend_z[npell*(nsd)*num_intb*k2+npell*(nsd)*k+npell*1+i] = 
		      *(dzdt_node+i)*(*(xs+nsd*2+1));
       		   shg.bend_z[npell*(nsd)*num_intb*k2+npell*(nsd)*k+npell*2+i] = 
		      *(dzdt_node+i)*(*(xs+nsd*2+2));

               	}
	   }
	}

/* 
   calculating the dN/dx,dy,dz and dza/dx,dy,dz matrices for 1X1 
   integration in lamina 
*/
/* 
   calculating dx,dy,dz/ds,dr matrix for 1X1 
   integration in lamina
*/
        for( j = 0; j < nsdl; ++j ) /* j loops over s,r */
	{
        	for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
        	{
	     	    check=dotX((xlxs.bar+nsdl*i+j),(shl.shear+npell*j),
			(xl.bar+npell*i), npell);
		}
	}
        for( k2 = 0; k2 < num_ints; ++k2 )
	{
               	for( j = 0; j < nsdl; ++j ) /* j loops over s,r */
	       	{
           	   for( i = 0; i < npell; ++i )
	   	   {
	      	   	*(shlshl2_vec+i)=
		   		shl.shear[npell*j+i]*(*(znode+npell*k2+i));
		   }

              	   for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
        	   {
	     	       check=dotX((xlxs.hat+nsd*i+j),shlshl2_vec,
				(xl.hat+npell*i), npell);
		       *(xs+nsd*i+j)=xlxs.bar[nsdl*i+j]+xlxs.hat[nsd*i+j];
		   }
	       	}
/* 
   calculating dx,dy,dz/dt matrix for 1X1 
   integration in lamina
*/
           	for( i = 0; i < npell; ++i )
	   	{
	      	   *(shlshl2_vec+i)=
			shl.shear[npell*2+i]*(*(dzdt_node+i));
		}
               	for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
               	{
	           check= dotX((xlxs.hat+nsd*i+2),shlshl2_vec,
			(xl.hat+npell*i), npell);
		   *(xs+nsd*i+2)=xlxs.hat[nsd*i+2];
	       	}

/* 
   calculating rotation matrix for lamina q[i,j] matrix for 1X1 
   integration in lamina
*/
		*(vecl_xi)=*(xs);		*(vecl_eta)=*(xs+1);
		*(vecl_xi+1)=*(xs+nsd*1);	*(vecl_eta+1)=*(xs+nsd*1+1);
		*(vecl_xi+2)=*(xs+nsd*2);	*(vecl_eta+2)=*(xs+nsd*2+1);

        	check = normal(vecl_xi);
		if(!check) printf( " Problems with vecl_xi element %d\n",el);
        	check = normal(vecl_eta);
		if(!check) printf( " Problems with vecl_eta element %d\n",el);
        	check = normcrossX(vecl_xi,vecl_eta,vecl_3);
		if(!check) printf( " Problems with normcrossX \n");

		*(vecl_alp)=*(vecl_xi)+*(vecl_eta);
		*(vecl_alp+1)=*(vecl_xi+1)+*(vecl_eta+1);
		*(vecl_alp+2)=*(vecl_xi+2)+*(vecl_eta+2);
        	check = normal(vecl_alp);
		if(!check) printf( " Problems with vecl_alp element %d\n",el);
        	check = normcrossX(vecl_3,vecl_alp,vecl_beta);
		if(!check) printf( " Problems with normcrossX \n");

		*(vecl_1)=rt2*(*(vecl_alp)-*(vecl_beta));
		*(vecl_1+1)=rt2*(*(vecl_alp+1)-*(vecl_beta+1));
		*(vecl_1+2)=rt2*(*(vecl_alp+2)-*(vecl_beta+2));
		*(vecl_2)=rt2*(*(vecl_alp)+*(vecl_beta));
		*(vecl_2+1)=rt2*(*(vecl_alp+1)+*(vecl_beta+1));
		*(vecl_2+2)=rt2*(*(vecl_alp+2)+*(vecl_beta+2));

               	for( i = 0; i < nsd; ++i ) 
               	{
		    rotate.l_shear[nsdsq*k2+i]=*(vecl_1+i);
		    rotate.l_shear[nsdsq*k2+nsd*1+i]=*(vecl_2+i);
		    rotate.l_shear[nsdsq*k2+nsd*2+i]=*(vecl_3+i);
		}

               	*(temp)=*(xs+4)*(*(xs+8))-*(xs+7)*(*(xs+5));
               	*(temp+3)=*(xs+6)*(*(xs+5))-*(xs+3)*(*(xs+8));
               	*(temp+6)=*(xs+3)*(*(xs+7))-*(xs+6)*(*(xs+4));
 
               	*(det+num_int+k2)=
			*(xs)*(*(temp))+*(xs+1)*(*(temp+3))+*(xs+2)*(*(temp+6));
               	/*printf("%d %f\n", num_intb, *(det+num_int+k2));*/
 
               	if(*(det+num_int+k2) <= 0.0 )
               	{
                    printf("the element (%d) is inverted; det:%f; fiber integ pt.:%d\n",
                        el,*(det+num_int+k2),k2);
                    return 0;
               	}
/* The inverse of the jacobian, dc/dx, is calculated below */
 
               	*(temp+1)=*(xs+7)*(*(xs+2))-*(xs+1)*(*(xs+8));
               	*(temp+4)=*(xs)*(*(xs+8))-*(xs+6)*(*(xs+2));
               	*(temp+7)=*(xs+6)*(*(xs+1))-*(xs)*(*(xs+7));
               	*(temp+2)=*(xs+1)*(*(xs+5))-*(xs+4)*(*(xs+2));
               	*(temp+5)=*(xs+3)*(*(xs+2))-*(xs)*(*(xs+5));
               	*(temp+8)=*(xs)*(*(xs+4))-*(xs+3)*(*(xs+1));
     
               	for( j = 0; j < nsd; ++j )
               	{
                    for( i = 0; i < nsd; ++i )
                    {
                       *(xs+nsd*i+j)=
			   *(temp+nsd*i+j)/(*(det+num_int+k2));
                    }
               	}
               	for( i = 0; i < npell; ++i )
               	{
                    *(col1)=shl.shear[i];
                    *(col1+1)=shl.shear[npell*1+i];
                    *(col2)=*(xs);
                    *(col2+1)=*(xs+nsd*1);
                    check=dotX((shg.shear+npell*(nsd+1)*k2+i),
			(col2),(col1),nsdl);
                    *(col2)=*(xs+1);
                    *(col2+1)=*(xs+nsd*1+1);
                    check=dotX((shg.shear+npell*(nsd+1)*k2+npell*1+i),
			(col2),(col1),nsdl);
                    *(col2)=*(xs+2);
                    *(col2+1)=*(xs+nsd*1+2);
                    check=dotX((shg.shear+npell*(nsd+1)*k2+npell*2+i),
			(col2),(col1),nsdl);
                    shg.shear[npell*(nsd+1)*k2+npell*3+i]=
			shl.shear[npell*2+i];

       		    shg.shear_z[npell*(nsd)*k2+i] =  
			*(dzdt_node+i)*(*(xs+nsd*2));
       		    shg.shear_z[npell*(nsd)*k2+npell*1+i] = 
			*(dzdt_node+i)*(*(xs+nsd*2+1));
       		    shg.shear_z[npell*(nsd)*k2+npell*2+i] = 
			*(dzdt_node+i)*(*(xs+nsd*2+2));
               	}
 
	}
        return 1; 
}


int shshg_mass( double *det, int el, SH shl, XL xl, double *zp1, double *zm1, 
	double *znode, double *dzdt_node)
{
/*
     This subroutine calculates the determinant for
     a shell element at the gauss points for the calculation of
     the global mass matrix.  Unlike shshg, we do not have to
     calculate the shape function derivatives with respect to
     the global coordinates x, y, z which are only needed for
     strains and stresses.

     It is based on the subroutine QDCSHG from the book
     "The Finite Element Method" by Thomas Hughes, page 783.
     This subroutine calculates the global shape function derivatives for
     a shell element at the gauss points.

 ....  CALCULATE GLOBAL DERIVATIVES OF SHAPE FUNCTIONS AND
       JACOBIAN DETERMINANTS FOR A FOUR-NODE QUADRALATERAL ELEMENT

       xl.bar[j+npell*i] = GLOBAL COORDINATES ON LAMINA 
       xl.hat[j+npell*i] = GLOBAL COORDINATES ON FIBER
       *(det+num_intb*k2+k)  = JACOBIAN DETERMINANT

       FOR 2X2 PT. GAUSS
       shl.bend[npell*(nsdl+1)*k+i] = LOCAL ("XI") DERIVATIVE OF SHAPE FUNCTION
       shl.bend[npell*(nsdl+1)*k+npell*1+i] = LOCAL ("ETA") DERIVATIVE OF SHAPE FUNCTION
       shl.bend[npell*(nsdl+1)*k+npell*2+i] = LOCAL SHAPE FUNCTION(in XI and ETA)

       FOR 2X2 PT. GAUSS
       shl.bend_z[i] = LOCAL ("ZETA") DERIVATIVE OF SHAPE FUNCTION
       shl.bend_z[npelf*num_ints+npelf*k2+i] = LOCAL SHAPE FUNCTION(in ZETA)

       *(xs+2*j+i) = JACOBIAN MATRIX
          i    = LOCAL NODE NUMBER OR GLOBAL COORDINATE NUMBER
          j    = GLOBAL COORDINATE NUMBER
          k    = INTEGRATION-POINT NUMBER FOR LAMINA 
          k2   = INTEGRATION-POINT NUMBER FOR FIBER
       num_intb    = NUMBER OF INTEGRATION POINTS FOR LAMINA, 4
       num_ints    = NUMBER OF INTEGRATION POINTS FOR FIBER, 2
	
			Updated 9/26/01
*/
        double xs[soxshat],temp[nsdsq],col1[nsdl],col2[nsdl],temp1,temp2;
	double shlshl2_vec[npell],xfib2[2];
	XLXS xlxs;
	int check,i,j,k,k2;

/* 
   calculating the dN/dx,dy,dz dza/dx,dy,dz matrices for 2X2 
   integration
*/
        for( k2 = 0; k2 < num_ints; ++k2 )
	{
           for( i = 0; i < npell; ++i )
	   {
/* 
   calculating the values of the za and dzdt matrix for 2X2 
   integration
*/
	       *(znode+npell*k2+i)=
		    shl.bend_z[npelf*(nsdf+1)*k2+npelf*1+1]*(*(zp1+i))+
		    shl.bend_z[npelf*(nsdf+1)*k2+npelf*1]*(*(zm1+i));
	       *(dzdt_node+i)=
		    shl.bend_z[npelf*(nsdf+1)*k2+1]*(*(zp1+i))+
		    shl.bend_z[npelf*(nsdf+1)*k2]*(*(zm1+i));
	   }
           for( k = 0; k < num_intb; ++k )
	   {
/* 
   calculating dx,dy,dz/ds,dr matrix for 2X2 
   integration
*/
           	for( j = 0; j < nsdl; ++j ) /* j loops over s,r */
	   	{
        		for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
        		{
	     		    check=dotX((xlxs.bar+nsdl*i+j),(shl.bend+npell*(nsdl+1)*k+npell*j),
				(xl.bar+npell*i),npell);
			}
	   	}
               	for( j = 0; j < nsdl; ++j ) /* j loops over s,r */
	       	{
           	   for( i = 0; i < npell; ++i )
	   	   {
	      	   	*(shlshl2_vec+i)=
			   shl.bend[npell*(nsdl+1)*k+npell*j+i]*(*(znode+npell*k2+i));
		   }
              	   for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
        	   {
	     	       check=dotX((xlxs.hat+nsd*i+j),shlshl2_vec,
				(xl.hat+npell*i), npell);
                       *(xs+nsd*i+j)=xlxs.bar[nsdl*i+j]+xlxs.hat[nsd*i+j];
		   }
	       	}
/* 
   calculating dx,dy,dz/dt matrix for 2X2 
   integration
*/
           	for( i = 0; i < npell; ++i )
	   	{
	      	   *(shlshl2_vec+i)=
			shl.bend[npell*(nsdl+1)*k+npell*2+i]*(*(dzdt_node+i));
		}
               	for( i = 0; i < nsd; ++i ) /* i loops over x,y,z */
               	{
	           check= dotX((xlxs.hat+nsd*i+2),shlshl2_vec,
			(xl.hat+npell*i), npell);
		   *(xs+nsd*i+2)=xlxs.hat[nsd*i+2];
	       	}

               	*(temp)=*(xs+4)*(*(xs+8))-*(xs+7)*(*(xs+5));
               	*(temp+3)=*(xs+6)*(*(xs+5))-*(xs+3)*(*(xs+8));
               	*(temp+6)=*(xs+3)*(*(xs+7))-*(xs+6)*(*(xs+4));

               	*(det+num_intb*k2+k)=
			*(xs)*(*(temp))+*(xs+1)*(*(temp+3))+*(xs+2)*(*(temp+6));
               	/*printf(" %d %f\n", k, *(det+num_intb*k2+k));*/
 
               	if(*(det+num_intb*k2+k) <= 0.0 )
               	{
                    printf("the element (%d) is inverted; det:%f; integ pt.:%d\n",
                        el,*(det+num_intb*k2+k),k);
                    return 0;
               	}
 
	   }
	}

        return 1; 
}

