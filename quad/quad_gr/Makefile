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


all: qdpost 

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_QD = ../quad/memory.o ../quad/qdreader.o 

OBJECTS_QDH = ../quad/qdconst.h ../quad/qdstruct.h 

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

PROGS = qdpost 

CC = /usr/bin/gcc

qdpost: qdpost.o \
		memgr.o \
		readgr2.o \
		param2.o \
		regetparam2.o \
		qdset.o \
		qdnorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		qdmesh.o \
		mshkey2.o \
		menu.o \
		qdgetnew.o \
		crtldsp2.o \
		crtlmse2.o \
		qdrender.o \
		$(OBJECTS_QDH) \
		$(OBJECTS_GRH) \
		qdgui.h \
		qdstrcgr.h 
	$(CC) -o qdpost $(LIB1) qdpost.o \
		memgr.o \
		readgr2.o \
		param2.o \
		regetparam2.o \
		qdset.o \
		qdnorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		qdmesh.o \
		mshkey2.o \
		menu.o \
		qdgetnew.o \
		crtldsp2.o \
		crtlmse2.o \
		qdrender.o \
		$(OBJECTS_QD) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

qdrender.o: qdrender.c 
	$(CC) -c qdrender.c 

crtlmse2.o: $(COMMON_GR)/crtlmse2.c 
	$(CC) -c $(COMMON_GR)/crtlmse2.c -DQUAD1

crtldsp2.o: $(COMMON_GR)/crtldsp2.c 
	$(CC) -c $(COMMON_GR)/crtldsp2.c -DQUAD1

qdgetnew.o: qdgetnew.c 
	$(CC) -c qdgetnew.c 

menu.o: $(COMMON_GR)/menu.c
	$(CC) -c $(COMMON_GR)/menu.c -DQUAD1

mshkey2.o: $(COMMON_GR)/mshkey2.c 
	$(CC) -c $(COMMON_GR)/mshkey2.c -DQUAD1

qdmesh.o: qdmesh.c 
	$(CC) -c qdmesh.c 

mshdsp.o: $(COMMON_GR)/mshdsp.c 
	$(CC) -c $(COMMON_GR)/mshdsp.c -DQUAD1

bc_vec.o: $(COMMON_GR)/bc_vec.c 
	$(CC) -c $(COMMON_GR)/bc_vec.c -DQUAD1

qdnorm_vec.o: qdnorm_vec.c 
	$(CC) -c qdnorm_vec.c -DQUAD1

qdset.o: qdset.c 
	$(CC) -c qdset.c 

regetparam2.o: $(COMMON_GR)/regetparam2.c 
	$(CC) -c $(COMMON_GR)/regetparam2.c -DQUAD1

param2.o: $(COMMON_GR)/param2.c 
	$(CC) -c $(COMMON_GR)/param2.c -DQUAD1

readgr2.o: $(COMMON_GR)/readgr2.c 
	$(CC) -c $(COMMON_GR)/readgr2.c -DQUAD1

memgr.o: $(COMMON_GR)/memgr.c 
	$(CC) -c $(COMMON_GR)/memgr.c -DQUAD1

qdpost.o: qdpost.c 
	$(CC) -c qdpost.c

clean:
	-rm $(PROGS)
	-rm *.o 

