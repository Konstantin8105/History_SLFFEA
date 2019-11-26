# SLFFEA Makefile
# Version:  1.5
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


all: br2post 

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_BR = ../brick2/memory.o ../brick2/br2reader.o 

OBJECTS_BRG = ../brick_gr/brnorm_vec.o 

OBJECTS_BRH = ../brick/brconst.h ../brick2/br2struct.h 

COMMON = ../../common

COMMON_GR = ../../common_gr

OBJECTS_GR = $(COMMON_GR)/pwinston.o $(COMMON_GR)/axes.o \
	$(COMMON_GR)/mshcommon.o $(COMMON_GR)/crtlcommon.o \
	$(COMMON_GR)/printext.o $(COMMON_GR)/screenshot.o \
	$(COMMON_GR)/filecheck.o $(COMMON_GR)/pointlocate.o \
	$(COMMON_GR)/axeslabel.o $(COMMON_GR)/quat.o \
	$(COMMON_GR)/crossplane.o 


OBJECTS_GRH = $(COMMON_GR)/color_gr.h $(COMMON_GR)/control.h 

OBJECTS2 = $(COMMON)/crossx.o $(COMMON)/matx.o 

PROGS = br2post 

CC = /usr/bin/gcc
LIBS = -lm

br2post: br2post.o \
		memgr.o \
		readgr.o \
		con_surf_read.o \
		param.o \
		regetparam.o \
		set.o \
		bc_vec.o \
		mshdsp.o \
		brmesh.o \
		mshkey.o \
		menu.o \
		br2getnew.o \
		crtldsp.o \
		crtlmse.o \
		brrender.o \
		$(OBJECTS_BRH) \
		$(OBJECTS_GRH) \
		br2gui.h 
	$(CC) -o br2post $(LIB1) br2post.o \
		memgr.o \
		readgr.o \
		con_surf_read.o \
		param.o \
		regetparam.o \
		set.o \
		bc_vec.o \
		mshdsp.o \
		brmesh.o \
		mshkey.o \
		menu.o \
		br2getnew.o \
		crtldsp.o \
		crtlmse.o \
		brrender.o \
		$(OBJECTS_BR) $(OBJECTS_BRG) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

brrender.o: ../brick_gr/brrender.c 
	$(CC) -c ../brick_gr/brrender.c -DBRICK2

crtlmse.o: $(COMMON_GR)/crtlmse.c 
	$(CC) -c $(COMMON_GR)/crtlmse.c -DBRICK2

crtldsp.o: $(COMMON_GR)/crtldsp.c 
	$(CC) -c $(COMMON_GR)/crtldsp.c -DBRICK2

br2getnew.o: br2getnew.c 
	$(CC) -c br2getnew.c 

menu.o: $(COMMON_GR)/menu.c 
	$(CC) -c $(COMMON_GR)/menu.c -DBRICK2

mshkey.o: $(COMMON_GR)/mshkey.c 
	$(CC) -c $(COMMON_GR)/mshkey.c -DBRICK2

brmesh.o: ../brick_gr/brmesh.c 
	$(CC) -c ../brick_gr/brmesh.c -DBRICK2

mshdsp.o: $(COMMON_GR)/mshdsp.c 
	$(CC) -c $(COMMON_GR)/mshdsp.c -DBRICK2

bc_vec.o: $(COMMON_GR)/bc_vec.c 
	$(CC) -c $(COMMON_GR)/bc_vec.c -DBRICK2

set.o: $(COMMON_GR)/set.c
	$(CC) -c $(COMMON_GR)/set.c -DBRICK2

regetparam.o: $(COMMON_GR)/regetparam.c
	$(CC) -c $(COMMON_GR)/regetparam.c -DBRICK2

param.o: $(COMMON_GR)/param.c
	$(CC) -c $(COMMON_GR)/param.c -DBRICK2

con_surf_read.o: $(COMMON_GR)/con_surf_read.c 
	$(CC) -c $(COMMON_GR)/con_surf_read.c -DBRICK2

readgr.o: $(COMMON_GR)/readgr.c 
	$(CC) -c $(COMMON_GR)/readgr.c -DBRICK2

memgr.o: $(COMMON_GR)/memgr.c 
	$(CC) -c $(COMMON_GR)/memgr.c -DBRICK2

br2post.o: br2post.c 
	$(CC) -c br2post.c

clean:
	-rm $(PROGS)
	-rm *.o 
