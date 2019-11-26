# SLFFEA Makefile
# Version:  1.5
# Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005  San Le (slffea@yahoo.com)
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


all: trpost 

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_QD = ../tri/memory2.o ../tri/trreader.o 

OBJECTS_QDH = ../tri/trconst.h ../tri/trstruct.h 

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

PROGS = trpost 

CC = /usr/bin/gcc

trpost: trpost.o \
		memgr2.o \
		readgr2.o \
		param2.o \
		regetparam2.o \
		trset.o \
		trnorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		trmesh.o \
		mshkey2.o \
		menu2.o \
		trgetnew.o \
		crtldsp2.o \
		crtlmse2.o \
		trrender.o \
		$(OBJECTS_QDH) \
		$(OBJECTS_GRH) \
		trgui.h \
		trstrcgr.h 
	$(CC) -o trpost $(LIB1) trpost.o \
		memgr2.o \
		readgr2.o \
		param2.o \
		regetparam2.o \
		trset.o \
		trnorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		trmesh.o \
		mshkey2.o \
		menu2.o \
		trgetnew.o \
		crtldsp2.o \
		crtlmse2.o \
		trrender.o \
		$(OBJECTS_QD) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

trrender.o: trrender.c 
	$(CC) -c trrender.c 

crtlmse2.o: $(COMMON_GR)/crtlmse2.c 
	$(CC) -c $(COMMON_GR)/crtlmse2.c -DTRI1

crtldsp2.o: $(COMMON_GR)/crtldsp2.c 
	$(CC) -c $(COMMON_GR)/crtldsp2.c -DTRI1

trgetnew.o: trgetnew.c 
	$(CC) -c trgetnew.c 

menu2.o: $(COMMON_GR)/menu2.c 
	$(CC) -c $(COMMON_GR)/menu2.c -DTRI1

mshkey2.o: $(COMMON_GR)/mshkey2.c 
	$(CC) -c $(COMMON_GR)/mshkey2.c -DTRI1

trmesh.o: trmesh.c 
	$(CC) -c trmesh.c 

mshdsp.o: $(COMMON_GR)/mshdsp.c 
	$(CC) -c $(COMMON_GR)/mshdsp.c -DTRI1

bc_vec.o: $(COMMON_GR)/bc_vec.c
	$(CC) -c $(COMMON_GR)/bc_vec.c -DTRI1

trnorm_vec.o: trnorm_vec.c
	$(CC) -c trnorm_vec.c -DTRI1

trset.o: trset.c 
	$(CC) -c trset.c 

regetparam2.o: $(COMMON_GR)/regetparam2.c 
	$(CC) -c $(COMMON_GR)/regetparam2.c -DTRI1

param2.o: $(COMMON_GR)/param2.c 
	$(CC) -c $(COMMON_GR)/param2.c -DTRI1

readgr2.o: $(COMMON_GR)/readgr2.c
	$(CC) -c $(COMMON_GR)/readgr2.c -DTRI1

memgr2.o: $(COMMON_GR)/memgr2.c 
	$(CC) -c $(COMMON_GR)/memgr2.c -DTRI1

trpost.o: trpost.c 
	$(CC) -c trpost.c

clean:
	-rm $(PROGS)
	-rm *.o 

