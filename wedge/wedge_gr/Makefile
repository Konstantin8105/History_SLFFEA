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


all: wepost

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_WE = ../wedge/memory.o ../wedge/wereader.o 

OBJECTS_WEH = ../wedge/weconst.h ../wedge/westruct.h 

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

PROGS = wepost

CC = /usr/bin/gcc -g
LIBS = -lm

wepost: wepost.o \
		memgr.o \
		readgr.o \
		con_surf_read.o \
		param.o \
		regetparam.o \
		set.o \
		wenorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		wemesh.o \
		mshkey.o \
		menu.o \
		wegetnew.o \
		crtldsp.o \
		crtlmse.o \
		werender.o \
		$(OBJECTS_WEH) \
		$(OBJECTS_GRH) \
		wegui.h \
		westrcgr.h 
	$(CC) -o wepost $(LIB1) wepost.o \
		memgr.o \
		readgr.o \
		con_surf_read.o \
		param.o \
		regetparam.o \
		set.o \
		wenorm_vec.o \
		bc_vec.o \
		mshdsp.o \
		wemesh.o \
		mshkey.o \
		menu.o \
		wegetnew.o \
		crtldsp.o \
		crtlmse.o \
		werender.o \
		$(OBJECTS_WE) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

werender.o: werender.c 
	$(CC) -c werender.c -DWEDGE1

crtlmse.o: $(COMMON_GR)/crtlmse.c 
	$(CC) -c $(COMMON_GR)/crtlmse.c -DWEDGE1

crtldsp.o: $(COMMON_GR)/crtldsp.c 
	$(CC) -c $(COMMON_GR)/crtldsp.c -DWEDGE1

wegetnew.o: wegetnew.c 
	$(CC) -c wegetnew.c 

menu.o: $(COMMON_GR)/menu.c 
	$(CC) -c $(COMMON_GR)/menu.c -DWEDGE1

mshkey.o: $(COMMON_GR)/mshkey.c 
	$(CC) -c $(COMMON_GR)/mshkey.c -DWEDGE1

wemesh.o: wemesh.c 
	$(CC) -c wemesh.c -DWEDGE1

mshdsp.o: $(COMMON_GR)/mshdsp.c 
	$(CC) -c $(COMMON_GR)/mshdsp.c -DWEDGE1

bc_vec.o: $(COMMON_GR)/bc_vec.c 
	$(CC) -c $(COMMON_GR)/bc_vec.c -DWEDGE1

wenorm_vec.o: wenorm_vec.c 
	$(CC) -c wenorm_vec.c -DWEDGE1

set.o: $(COMMON_GR)/set.c 
	$(CC) -c $(COMMON_GR)/set.c -DWEDGE1

regetparam.o: $(COMMON_GR)/regetparam.c 
	$(CC) -c $(COMMON_GR)/regetparam.c -DWEDGE1

param.o: $(COMMON_GR)/param.c 
	$(CC) -c $(COMMON_GR)/param.c -DWEDGE1

con_surf_read.o: $(COMMON_GR)/con_surf_read.c
	$(CC) -c $(COMMON_GR)/con_surf_read.c -DWEDGE1

readgr.o: $(COMMON_GR)/readgr.c
	$(CC) -c $(COMMON_GR)/readgr.c -DWEDGE1

memgr.o: $(COMMON_GR)/memgr.c 
	$(CC) -c $(COMMON_GR)/memgr.c -DWEDGE1

wepost.o: wepost.c 
	$(CC) -c wepost.c

clean:
	-rm $(PROGS)
	-rm *.o 
