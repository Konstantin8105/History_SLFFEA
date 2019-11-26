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


all: trpost 

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_QD = ../tri/trmemory.o ../tri/trreader.o 

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
		trmemgr.o \
		trreadgr.o \
		trparam.o \
		trparam2.o \
		trset.o \
		trbc_vec.o \
		trmshdsp.o \
		trmesh.o \
		trmshkey.o \
		trmenu.o \
		trgetnew.o \
		trcrtldsp.o \
		trcrtlmse.o \
		trrender.o \
		$(OBJECTS_QDH) \
		$(OBJECTS_GRH) \
		trgui.h \
		trstrcgr.h 
	$(CC) -o trpost $(LIB1) trpost.o \
		trmemgr.o \
		trreadgr.o \
		trparam.o \
		trparam2.o \
		trset.o \
		trbc_vec.o \
		trmshdsp.o \
		trmesh.o \
		trmshkey.o \
		trmenu.o \
		trgetnew.o \
		trcrtldsp.o \
		trcrtlmse.o \
		trrender.o \
		$(OBJECTS_QD) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

trrender.o: trrender.c 
	$(CC) -c trrender.c 

trcrtlmse.o: trcrtlmse.c 
	$(CC) -c trcrtlmse.c 

trcrtldsp.o: trcrtldsp.c 
	$(CC) -c trcrtldsp.c 

trgetnew.o: trgetnew.c 
	$(CC) -c trgetnew.c 

trmenu.o: trmenu.c 
	$(CC) -c trmenu.c 

trmshkey.o: trmshkey.c 
	$(CC) -c trmshkey.c 

trmesh.o: trmesh.c 
	$(CC) -c trmesh.c 

trmshdsp.o: trmshdsp.c 
	$(CC) -c trmshdsp.c 

trbc_vec.o: trbc_vec.c 
	$(CC) -c trbc_vec.c 

trset.o: trset.c 
	$(CC) -c trset.c 

trparam2.o: trparam2.c 
	$(CC) -c trparam2.c 

trparam.o: trparam.c 
	$(CC) -c trparam.c 

trreadgr.o: trreadgr.c 
	$(CC) -c trreadgr.c 

trmemgr.o: trmemgr.c 
	$(CC) -c trmemgr.c 

trpost.o: trpost.c 
	$(CC) -c trpost.c

clean:
	-rm $(PROGS)
	-rm *.o 

