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


all: br2post 

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_BR = ../brick2/brmemory.o ../brick2/br2reader.o 

OBJECTS_BRG = ../brick_gr/brnorm_vec.o 

OBJECTS_BRH = ../brick/brconst.h ../brick2/br2struct.h 

COMMON = ../../common

COMMON_GR = ../../common_gr

OBJECTS_GR = $(COMMON_GR)/pwinston.o $(COMMON_GR)/axes.o \
	$(COMMON_GR)/mshcommon.o $(COMMON_GR)/crtlcommon.o \
	$(COMMON_GR)/printext.o $(COMMON_GR)/screenshot.o \
	$(COMMON_GR)/filecheck.o 

OBJECTS_GRH = $(COMMON_GR)/color_gr.h $(COMMON_GR)/control.h 

OBJECTS2 = $(COMMON)/crossx.o 

PROGS = br2post 

CC = /usr/bin/gcc
LIBS = -lm

br2post: br2post.o \
		brmemgr.o \
		brreadgr.o \
		brreadgr2.o \
		br2param.o \
		br2param2.o \
		br2set.o \
		brbc_vec.o \
		brmshdsp.o \
		brmesh.o \
		br2mshkey.o \
		br2menu.o \
		br2getnew.o \
		brcrtldsp.o \
		br2crtlmse.o \
		brrender.o \
		$(OBJECTS_BRH) \
		$(OBJECTS_GRH) \
		br2gui.h 
	$(CC) -o br2post $(LIB1) br2post.o \
		brmemgr.o \
		brreadgr.o \
		brreadgr2.o \
		br2param.o \
		br2param2.o \
		br2set.o \
		brbc_vec.o \
		brmshdsp.o \
		brmesh.o \
		br2mshkey.o \
		br2menu.o \
		br2getnew.o \
		brcrtldsp.o \
		br2crtlmse.o \
		brrender.o \
		$(OBJECTS_BR) $(OBJECTS_BRG) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

brrender.o: ../brick_gr/brrender.c 
	$(CC) -c ../brick_gr/brrender.c -DBRICK2

br2crtlmse.o: br2crtlmse.c 
	$(CC) -c br2crtlmse.c 

brcrtldsp.o: ../brick_gr/brcrtldsp.c 
	$(CC) -c ../brick_gr/brcrtldsp.c -DBRICK2

br2getnew.o: br2getnew.c 
	$(CC) -c br2getnew.c 

br2menu.o: br2menu.c 
	$(CC) -c br2menu.c 

br2mshkey.o: br2mshkey.c 
	$(CC) -c br2mshkey.c 

brmesh.o: ../brick_gr/brmesh.c 
	$(CC) -c ../brick_gr/brmesh.c -DBRICK2

brmshdsp.o: ../brick_gr/brmshdsp.c 
	$(CC) -c ../brick_gr/brmshdsp.c -DBRICK2

brbc_vec.o: ../brick_gr/brbc_vec.c 
	$(CC) -c ../brick_gr/brbc_vec.c -DBRICK2

br2set.o: br2set.c 
	$(CC) -c br2set.c 

br2param2.o: br2param2.c 
	$(CC) -c br2param2.c 

br2param.o: br2param.c 
	$(CC) -c br2param.c 

brreadgr2.o: ../brick_gr/brreadgr2.c 
	$(CC) -c ../brick_gr/brreadgr2.c -DBRICK2

brreadgr.o: ../brick_gr/brreadgr.c 
	$(CC) -c ../brick_gr/brreadgr.c -DBRICK2

brmemgr.o: ../brick_gr/brmemgr.c 
	$(CC) -c ../brick_gr/brmemgr.c -DBRICK2

br2post.o: br2post.c 
	$(CC) -c br2post.c

clean:
	-rm $(PROGS)
	-rm *.o 
