# SLFFEA Makefile
# Version:  1.2
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

#May work on Sun
#C = /usr/lang/gcc  -g

#Dubugging Compile for Linux
CC = /usr/bin/gcc  -g

#Optimization Compile for Linux
#CC = /usr/bin/gcc  -O3

LIBS = -lm
COMMON = ../../common
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/dotx.o \
	$(COMMON)/lm.o $(COMMON)/idiag.o $(COMMON)/kasmbl.o \
	$(COMMON)/sky.o  $(COMMON)/qr.o $(COMMON)/eigen.o \
	$(COMMON)/eigenprint.o $(COMMON)/rotate.o

PROGS = bm bmp

all: bm bmp

bm: fembeam.o \
		bmmemory.o \
		bmreader.o \
		bmsort.o \
		bmid.o \
		bmshape.o \
		bmcrossx.o \
		bmkasmbl.o \
		bmconj.o \
		bmbound.o \
		bmmasmbl.o \
		bmlanczos.o \
		bmmassp.o \
		bmwriter.o \
		bmshape_struct.h \
		bmconst.h \
		bmstruct.h 
	$(CC) -o bm fembeam.o \
		bmmemory.o \
		bmreader.o \
		bmsort.o \
		bmid.o \
		bmshape.o \
		bmcrossx.o \
		bmkasmbl.o \
		bmconj.o \
		bmbound.o \
		bmmasmbl.o \
		bmlanczos.o \
		bmmassp.o \
		bmwriter.o $(OBJECTS) $(LIBS)

bmp: bmp.o \
		bmpatch.o
	$(CC) -o bmp bmpatch.o $(LIBS)

bmp.o: bmpatch.c
	$(CC) -c bmpatch.c

bmwriter.o: bmwriter.c
	$(CC) -c bmwriter.c 

bmmassp.o: bmmassp.c
	$(CC) -c bmmassp.c 

bmlanczos.o: bmlanczos.c
	$(CC) -c bmlanczos.c 

bmmasmbl.o: bmmasmbl.c
	$(CC) -c bmmasmbl.c 

bmbound.o: bmbound.c
	$(CC) -c bmbound.c 

bmconj.o: bmconj.c
	$(CC) -c bmconj.c 

bmkasmbl.o: bmkasmbl.c
	$(CC) -c bmkasmbl.c 

bmcrossx.o: bmcrossx.c
	$(CC) -c bmcrossx.c 

bmshape.o: bmshape.c
	$(CC) -c bmshape.c 

bmid.o: bmid.c
	$(CC) -c bmid.c 

bmsort.o: bmsort.c
	$(CC) -c bmsort.c

bmreader.o: bmreader.c
	$(CC) -c bmreader.c

bmmemory.o: bmmemory.c
	$(CC) -c bmmemory.c

fembeam.o: fembeam.c
	$(CC) -c fembeam.c -DB3D

clean:
	-rm *.o $(PROGS)
