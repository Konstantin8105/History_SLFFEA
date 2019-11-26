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


all: bmpost 

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_BM = ../beam/bmmemory.o ../beam/bmreader.o 

OBJECTS_BMH = ../beam/bmconst.h ../beam/bmstruct.h 

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

COMMON_GR = ../../common_gr

PROGS = bmpost 

CC = /usr/bin/gcc -g

bmpost: bmpost.o \
		bmmemgr.o \
		bmreadgr.o \
		bmparam.o \
		bmregetparam.o \
		bmset.o \
		bmbc_vec.o \
		bmrotate.o \
		bmmshdsp.o \
		bmmesh.o \
		bmmshkey.o \
		bmmshkey2.o \
		bmmenu.o \
		bmgetnew.o \
		bmcrtldsp.o \
		bmcrtlmse.o \
		$(OBJECTS_BMH) \
		$(OBJECTS_GRH) \
		bmgui.h \
		bmstrcgr.h 
	$(CC) -o bmpost $(LIB1) bmpost.o \
		bmmemgr.o \
		bmreadgr.o \
		bmparam.o \
		bmregetparam.o \
		bmset.o \
		bmbc_vec.o \
		bmrotate.o \
		bmmshdsp.o \
		bmmesh.o \
		bmmshkey.o \
		bmmshkey2.o \
		bmmenu.o \
		bmgetnew.o \
		bmcrtldsp.o \
		bmcrtlmse.o \
		$(OBJECTS_BM) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

bmcrtlmse.o: bmcrtlmse.c 
	$(CC) -c bmcrtlmse.c 

bmcrtldsp.o: bmcrtldsp.c 
	$(CC) -c bmcrtldsp.c 

bmgetnew.o: bmgetnew.c 
	$(CC) -c bmgetnew.c 

bmmenu.o: bmmenu.c 
	$(CC) -c bmmenu.c 

bmmshkey2.o: bmmshkey2.c 
	$(CC) -c bmmshkey2.c 

bmmshkey.o: bmmshkey.c 
	$(CC) -c bmmshkey.c 

bmmesh.o: bmmesh.c 
	$(CC) -c bmmesh.c 

bmmshdsp.o: bmmshdsp.c 
	$(CC) -c bmmshdsp.c 

bmrotate.o: bmrotate.c 
	$(CC) -c bmrotate.c 

bmbc_vec.o: bmbc_vec.c 
	$(CC) -c bmbc_vec.c 

bmset.o: bmset.c 
	$(CC) -c bmset.c 

bmregetparam.o: bmregetparam.c 
	$(CC) -c bmregetparam.c 

bmparam.o: bmparam.c 
	$(CC) -c bmparam.c 

bmreadgr.o: bmreadgr.c 
	$(CC) -c bmreadgr.c 

bmmemgr.o: bmmemgr.c 
	$(CC) -c bmmemgr.c 

bmpost.o: bmpost.c 
	$(CC) -c bmpost.c

clean:
	-rm $(PROGS)
	-rm *.o 

targets: $(PROGS)

