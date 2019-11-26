# SLFFEA Makefile
# Version:  1.3
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

OBJECTS_QD = ../quad/qdmemory.o ../quad/qdreader.o 

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
		qdmemgr.o \
		qdreadgr.o \
		qdparam.o \
		qdparam2.o \
		qdset.o \
		qdbc_vec.o \
		qdmshdsp.o \
		qdmesh.o \
		qdmshkey.o \
		qdmenu.o \
		qdgetnew.o \
		qdcrtldsp.o \
		qdcrtlmse.o \
		qdrender.o \
		$(OBJECTS_QDH) \
		$(OBJECTS_GRH) \
		qdgui.h \
		qdstrcgr.h 
	$(CC) -o qdpost $(LIB1) qdpost.o \
		qdmemgr.o \
		qdreadgr.o \
		qdparam.o \
		qdparam2.o \
		qdset.o \
		qdbc_vec.o \
		qdmshdsp.o \
		qdmesh.o \
		qdmshkey.o \
		qdmenu.o \
		qdgetnew.o \
		qdcrtldsp.o \
		qdcrtlmse.o \
		qdrender.o \
		$(OBJECTS_QD) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

qdrender.o: qdrender.c 
	$(CC) -c qdrender.c 

qdcrtlmse.o: qdcrtlmse.c 
	$(CC) -c qdcrtlmse.c 

qdcrtldsp.o: qdcrtldsp.c 
	$(CC) -c qdcrtldsp.c 

qdgetnew.o: qdgetnew.c 
	$(CC) -c qdgetnew.c 

qdmenu.o: qdmenu.c 
	$(CC) -c qdmenu.c 

qdmshkey.o: qdmshkey.c 
	$(CC) -c qdmshkey.c 

qdmesh.o: qdmesh.c 
	$(CC) -c qdmesh.c 

qdmshdsp.o: qdmshdsp.c 
	$(CC) -c qdmshdsp.c 

qdbc_vec.o: qdbc_vec.c 
	$(CC) -c qdbc_vec.c 

qdset.o: qdset.c 
	$(CC) -c qdset.c 

qdparam2.o: qdparam2.c 
	$(CC) -c qdparam2.c 

qdparam.o: qdparam.c 
	$(CC) -c qdparam.c 

qdreadgr.o: qdreadgr.c 
	$(CC) -c qdreadgr.c 

qdmemgr.o: qdmemgr.c 
	$(CC) -c qdmemgr.c 

qdpost.o: qdpost.c 
	$(CC) -c qdpost.c

clean:
	-rm $(PROGS)
	-rm *.o 

