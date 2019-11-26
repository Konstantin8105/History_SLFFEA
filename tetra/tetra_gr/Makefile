# SLFFEA Makefile
# Version:  1.4
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


all: tepost

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_TE = ../tetra/memory2.o ../tetra/tereader.o 

OBJECTS_TEH = ../tetra/teconst.h ../tetra/testruct.h 

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

PROGS = tepost

CC = /usr/bin/gcc
LIBS = -lm

tepost: tepost.o \
		memgr2.o \
		readgr.o \
		con_surf_read.o \
		param.o \
		regetparam.o \
		teset.o \
		tenorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		temesh.o \
		mshkey.o \
		menu2.o \
		tegetnew.o \
		crtldsp.o \
		crtlmse.o \
		terender.o \
		$(OBJECTS_TEH) \
		$(OBJECTS_GRH) \
		tegui.h \
		testrcgr.h 
	$(CC) -o tepost $(LIB1) tepost.o \
		memgr2.o \
		readgr.o \
		con_surf_read.o \
		param.o \
		regetparam.o \
		teset.o \
		tenorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		temesh.o \
		mshkey.o \
		menu2.o \
		tegetnew.o \
		crtldsp.o \
		crtlmse.o \
		terender.o \
		$(OBJECTS_TE) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

terender.o: terender.c 
	$(CC) -c terender.c -DTETRA1

crtlmse.o: $(COMMON_GR)/crtlmse.c 
	$(CC) -c $(COMMON_GR)/crtlmse.c -DTETRA1

crtldsp.o: $(COMMON_GR)/crtldsp.c 
	$(CC) -c $(COMMON_GR)/crtldsp.c -DTETRA1

tegetnew.o: tegetnew.c 
	$(CC) -c tegetnew.c 

menu2.o: $(COMMON_GR)/menu2.c 
	$(CC) -c $(COMMON_GR)/menu2.c -DTETRA1

mshkey.o: $(COMMON_GR)/mshkey.c 
	$(CC) -c $(COMMON_GR)/mshkey.c -DTETRA1

temesh.o: temesh.c 
	$(CC) -c temesh.c -DTETRA1

mshdsp.o: $(COMMON_GR)/mshdsp.c 
	$(CC) -c $(COMMON_GR)/mshdsp.c -DTETRA1

bc_vec.o: $(COMMON_GR)/bc_vec.c 
	$(CC) -c $(COMMON_GR)/bc_vec.c -DTETRA1

tenorm_vec.o: tenorm_vec.c 
	$(CC) -c tenorm_vec.c 

teset.o: teset.c 
	$(CC) -c teset.c 

regetparam.o: $(COMMON_GR)/regetparam.c 
	$(CC) -c $(COMMON_GR)/regetparam.c -DTETRA1

param.o: $(COMMON_GR)/param.c 
	$(CC) -c $(COMMON_GR)/param.c -DTETRA1

con_surf_read.o: $(COMMON_GR)/con_surf_read.c
	$(CC) -c $(COMMON_GR)/con_surf_read.c -DTETRA1

readgr.o: $(COMMON_GR)/readgr.c 
	$(CC) -c $(COMMON_GR)/readgr.c -DTETRA1

memgr2.o: $(COMMON_GR)/memgr2.c 
	$(CC) -c $(COMMON_GR)/memgr2.c -DTETRA1

tepost.o: tepost.c 
	$(CC) -c tepost.c

clean:
	-rm $(PROGS)
	-rm *.o 
