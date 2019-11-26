/*
    This program does a screen dump of the Mesh Window.  It
    is almost entirely based on a Usenet posting of code
    by Antony Searle on the group:
  
        comp.graphics.api.opengl
  
                                  San Le
  
                        Last Update 5/14/00
  
    You can reach him at:
  
    Antony Searle
    H1-NF National Plasma Fusion Research Facility
    Australian National University
    acs654@my-dejanews.com
  
    SLFFEA source file
    Version:  1.1

    The source code contained in this file is released under the
    terms of the GNU Library General Public License.
 
*/

#if WINDOWS
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void ScreenShot( int width, int height)
{
   unsigned char *ScreenBuffer;
    FILE *Handle;
    unsigned char Header[18];
	int dum;

/* The width and height have to be multiples of 20.  Targa seems to
   require this, or else rescaling of the box will result in
   distortion of the image
*/
	dum = width%4;
	width -=  dum;
	dum = height%4;
	height -= dum;

        /* use glReadPixels and save a .tga file */
        ScreenBuffer = (unsigned char *)
                calloc(3*4*width*height,sizeof(unsigned char));
        glReadPixels(0, 0, width, height, GL_BGR,
                GL_UNSIGNED_BYTE, ScreenBuffer);

    Header[ 0] = 0;
    Header[ 1] = 0;
    Header[ 2] = 2;     /* Uncompressed, uninteresting */
    Header[ 3] = 0;
    Header[ 4] = 0;
    Header[ 5] = 0;
    Header[ 6] = 0;
    Header[ 7] = 0;
    Header[ 8] = 0;
    Header[ 9] = 0;
    Header[10] = 0;
    Header[11] = 0;
    Header[12] = (unsigned char) width;  /* Dimensions */
    Header[13] = (unsigned char) ((unsigned long) width >> 8);
    Header[14] = (unsigned char) height;
    Header[15] = (unsigned char) ((unsigned long) height >> 8);
    Header[16] = 24;    /* Bits per pixel */
    Header[17] = 0;

    Handle = fopen("Screen.tga", "wb");
    if(Handle == NULL) {
                free(ScreenBuffer);
        return;
        }

    fseek(Handle, 0, 0);
    fwrite(Header, 1, 18, Handle);
    fseek(Handle, 18, 0);
    fwrite(ScreenBuffer, 3, width * height, Handle);
    fclose(Handle);

    free(ScreenBuffer);
}
