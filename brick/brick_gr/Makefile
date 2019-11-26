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


all: brpost pvbr

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_BR = ../brick/brmemory.o ../brick/brreader.o 

OBJECTS_BRH = ../brick/brconst.h ../brick/brstruct.h 

COMMON = ../../common

COMMON_GR = ../../common_gr

OBJECTS_GR = $(COMMON_GR)/pwinston.o $(COMMON_GR)/axes.o \
	$(COMMON_GR)/mshcommon.o $(COMMON_GR)/crtlcommon.o \
	$(COMMON_GR)/printext.o $(COMMON_GR)/screenshot.o \
	$(COMMON_GR)/filecheck.o

OBJECTS_GRH = $(COMMON_GR)/color_gr.h $(COMMON_GR)/control.h 

OBJECTS2 = $(COMMON)/crossx.o 

PROGS = brpost pvbr

CC = /usr/bin/gcc
LIBS = -lm

brpost: brpost.o \
		brmemgr.o \
		brreadgr.o \
		brreadgr2.o \
		brparam.o \
		brparam2.o \
		brset.o \
		brnorm_vec.o \
		brbc_vec.o \
		brmshdsp.o \
		brmesh.o \
		brmshkey.o \
		brmenu.o \
		brgetnew.o \
		brcrtldsp.o \
		brcrtlmse.o \
		brrender.o \
		$(OBJECTS_BRH) \
		$(OBJECTS_GRH) \
		brgui.h \
		brstrcgr.h 
	$(CC) -o brpost $(LIB1) brpost.o \
		brmemgr.o \
		brreadgr.o \
		brreadgr2.o \
		brparam.o \
		brparam2.o \
		brset.o \
		brnorm_vec.o \
		brbc_vec.o \
		brmshdsp.o \
		brmesh.o \
		brmshkey.o \
		brmenu.o \
		brgetnew.o \
		brcrtldsp.o \
		brcrtlmse.o \
		brrender.o \
		$(OBJECTS_BR) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

pvbr: pvbr.o \
		$(OBJECTS_BRH) 
	$(CC) -o pvbr brinvent.o \
		$(OBJECTS_BR)  $(OBJECTS2) $(LIBS)

pvbr.o: brinvent.c
	$(CC) -c brinvent.c

brrender.o: brrender.c 
	$(CC) -c brrender.c -DBRICK1

brcrtlmse.o: brcrtlmse.c 
	$(CC) -c brcrtlmse.c 

brcrtldsp.o: brcrtldsp.c 
	$(CC) -c brcrtldsp.c -DBRICK1

brgetnew.o: brgetnew.c 
	$(CC) -c brgetnew.c 

brmenu.o: brmenu.c 
	$(CC) -c brmenu.c 

brmshkey.o: brmshkey.c 
	$(CC) -c brmshkey.c 

brmesh.o: brmesh.c 
	$(CC) -c brmesh.c -DBRICK1

brmshdsp.o: brmshdsp.c 
	$(CC) -c brmshdsp.c -DBRICK1

brbc_vec.o: brbc_vec.c 
	$(CC) -c brbc_vec.c -DBRICK1

brnorm_vec.o: brnorm_vec.c 
	$(CC) -c brnorm_vec.c 

brset.o: brset.c 
	$(CC) -c brset.c 

brparam2.o: brparam2.c 
	$(CC) -c brparam2.c 

brparam.o: brparam.c 
	$(CC) -c brparam.c 

brreadgr2.o: brreadgr2.c 
	$(CC) -c brreadgr2.c -DBRICK1

brreadgr.o: brreadgr.c 
	$(CC) -c brreadgr.c -DBRICK1

brmemgr.o: brmemgr.c 
	$(CC) -c brmemgr.c -DBRICK1

brpost.o: brpost.c 
	$(CC) -c brpost.c

clean:
	-rm $(PROGS)
	-rm *.o 
