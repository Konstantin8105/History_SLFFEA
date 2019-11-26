# SLFFEA Makefile
# Version:  1.1
# Copyright (C) 1999, 2000  San Le (slffea@yahoo.com)
# http://slffea.com/
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the Free
# Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.



all: shpost pvs

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_SH = ../shell/shmemory.o ../shell/shreader.o 

OBJECTS_SHH = ../shell/shconst.h ../shell/shstruct.h 

COMMON = ../../common

COMMON_GR = ../../common_gr

OBJECTS_GR = $(COMMON_GR)/pwinston.o $(COMMON_GR)/axes.o \
	$(COMMON_GR)/mshcommon.o $(COMMON_GR)/crtlcommon.o \
	$(COMMON_GR)/printext.o $(COMMON_GR)/screenshot.o \
	$(COMMON_GR)/filecheck.o 

OBJECTS_GRH = $(COMMON_GR)/color_gr.h $(COMMON_GR)/control.h 

OBJECTS2 = $(COMMON)/crossx.o $(COMMON)/dotx.o 

COMMON_GR = ../../common_gr

PROGS = shpost pvs

CC = /usr/bin/gcc

LIBS = -lm

shpost: shpost.o \
		shmemgr.o \
		shreadgr.o \
		shparam.o \
		shparam2.o \
		shset.o \
		shnorm_vec.o \
		shbc_vec.o \
		shmshdsp.o \
		shmesh.o \
		shmshkey.o \
		shmenu.o \
		shgetnew.o \
		shcrtldsp.o \
		shcrtlmse.o \
		shrender.o \
		$(OBJECTS_SHH) \
		$(OBJECTS_GRH) \
		shgui.h \
		shstrcgr.h 
	$(CC) -o shpost $(LIB1) shpost.o \
		shmemgr.o \
		shreadgr.o \
		shparam.o \
		shparam2.o \
		shset.o \
		shnorm_vec.o \
		shbc_vec.o \
		shmshdsp.o \
		shmesh.o \
		shmshkey.o \
		shmenu.o \
		shgetnew.o \
		shcrtldsp.o \
		shcrtlmse.o \
		shrender.o \
		$(OBJECTS_SH) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

pvs: pvs.o \
		$(OBJECTS_SHH)
	$(CC) -o pvs shinvent.o \
		$(OBJECTS_SH)  $(OBJECTS2) $(LIBS)

pvs.o: shinvent.c
	$(CC) -c shinvent.c

shrender.o: shrender.c 
	$(CC) -c shrender.c 

shcrtlmse.o: shcrtlmse.c 
	$(CC) -c shcrtlmse.c 

shcrtldsp.o: shcrtldsp.c 
	$(CC) -c shcrtldsp.c 

shgetnew.o: shgetnew.c 
	$(CC) -c shgetnew.c 

shmenu.o: shmenu.c 
	$(CC) -c shmenu.c 

shmshkey.o: shmshkey.c 
	$(CC) -c shmshkey.c 

shmesh.o: shmesh.c 
	$(CC) -c shmesh.c 

shmshdsp.o: shmshdsp.c 
	$(CC) -c shmshdsp.c 

shbc_vec.o: shbc_vec.c 
	$(CC) -c shbc_vec.c 

shnorm_vec.o: shnorm_vec.c 
	$(CC) -c shnorm_vec.c 

shset.o: shset.c 
	$(CC) -c shset.c

shparam2.o: shparam2.c 
	$(CC) -c shparam2.c 

shparam.o: shparam.c 
	$(CC) -c shparam.c 

shreadgr.o: shreadgr.c 
	$(CC) -c shreadgr.c 

shmemgr.o: shmemgr.c 
	$(CC) -c shmemgr.c 

shpost.o: shpost.c 
	$(CC) -c shpost.c

clean:
	-rm $(PROGS)
	-rm *.o 

