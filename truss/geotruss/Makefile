# SLFFEA Makefile
# Version:  1.0
# Copyright (C) 1999  San Le (slffea.juno.com)
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

#May work on Sun
#C = /usr/lang/gcc  -g

#Dubugging Compile for Linux
CC = /usr/bin/gcc  -g

#Optimization Compile for Linux
#CC = /usr/bin/gcc  -O3

LIBS = -lm
PROGS = ft

all: ft 

ft: femtruss.o \
		tsmemory.o \
		tsreader.o \
		tsidlm.o \
		tsidiag.o \
		tskasmbl.o \
		matx.o \
		matxt.o \
		sky.o \
		dotx.o \
		tswriter.o \
		tsconst.h \
		tsstruct.h 
	$(CC) -o ft femtruss.o \
		tsmemory.o \
		tsreader.o \
		tsidlm.o \
		tsidiag.o \
		tskasmbl.o \
		matx.o \
		matxt.o \
		sky.o \
		dotx.o \
		tswriter.o $(LIBS)

tswriter.o: tswriter.c
	$(CC) -c tswriter.c 

matx.o: matx.c
	$(CC) -c matx.c 

matxt.o: matxt.c
	$(CC) -c matxt.c 

dotx.o: dotx.c
	$(CC) -c dotx.c 

sky.o: sky.c
	$(CC) -c sky.c 

femtruss.o: femtruss.c
	$(CC) -o femtruss.o -c femtruss.c -DDATA_ON

tskasmbl.o: tskasmbl.c
	$(CC) -o tskasmbl.o -c tskasmbl.c -DDATA_ON

tsreader.o: tsreader.c
	$(CC) -c tsreader.c

tsidiag.o: tsidiag.c
	$(CC) -c tsidiag.c 

tsidlm.o: tsidlm.c
	$(CC) -c tsidlm.c 

tsmemory.o: tsmemory.c
	$(CC) -c tsmemory.c

clean:
	-rm *.o $(PROGS)
