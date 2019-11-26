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


all: tepost

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_TE = ../tetra/tememory.o ../tetra/tereader.o 

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
		tememgr.o \
		tereadgr.o \
		tereadgr2.o \
		teparam.o \
		teparam2.o \
		teset.o \
		tenorm_vec.o \
		tebc_vec.o \
		temshdsp.o \
		temesh.o \
		temshkey.o \
		temenu.o \
		tegetnew.o \
		tecrtldsp.o \
		tecrtlmse.o \
		terender.o \
		$(OBJECTS_TEH) \
		$(OBJECTS_GRH) \
		tegui.h \
		testrcgr.h 
	$(CC) -o tepost $(LIB1) tepost.o \
		tememgr.o \
		tereadgr.o \
		tereadgr2.o \
		teparam.o \
		teparam2.o \
		teset.o \
		tenorm_vec.o \
		tebc_vec.o \
		temshdsp.o \
		temesh.o \
		temshkey.o \
		temenu.o \
		tegetnew.o \
		tecrtldsp.o \
		tecrtlmse.o \
		terender.o \
		$(OBJECTS_TE) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

terender.o: terender.c 
	$(CC) -c terender.c -DTETRA1

tecrtlmse.o: tecrtlmse.c 
	$(CC) -c tecrtlmse.c 

tecrtldsp.o: tecrtldsp.c 
	$(CC) -c tecrtldsp.c -DTETRA1

tegetnew.o: tegetnew.c 
	$(CC) -c tegetnew.c 

temenu.o: temenu.c 
	$(CC) -c temenu.c 

temshkey.o: temshkey.c 
	$(CC) -c temshkey.c 

temesh.o: temesh.c 
	$(CC) -c temesh.c -DTETRA1

temshdsp.o: temshdsp.c 
	$(CC) -c temshdsp.c -DTETRA1

tebc_vec.o: tebc_vec.c 
	$(CC) -c tebc_vec.c -DTETRA1

tenorm_vec.o: tenorm_vec.c 
	$(CC) -c tenorm_vec.c 

teset.o: teset.c 
	$(CC) -c teset.c 

teparam2.o: teparam2.c 
	$(CC) -c teparam2.c 

teparam.o: teparam.c 
	$(CC) -c teparam.c 

tereadgr2.o: tereadgr2.c 
	$(CC) -c tereadgr2.c -DTETRA1

tereadgr.o: tereadgr.c 
	$(CC) -c tereadgr.c -DTETRA1

tememgr.o: tememgr.c 
	$(CC) -c tememgr.c -DTETRA1

tepost.o: tepost.c 
	$(CC) -c tepost.c

clean:
	-rm $(PROGS)
	-rm *.o 
