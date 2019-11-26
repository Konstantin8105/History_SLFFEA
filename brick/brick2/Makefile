# SLFFEA Makefile
# Version:  1.0
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
#CC = /usr/bin/gcc  -g

#Optimization Compile for Linux
CC = /usr/bin/gcc  -O3

LIBS = -lm
COMMON = ../../common
BRICK = ../brick
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/dotx.o \
	$(COMMON)/lm.o $(COMMON)/idiag.o $(COMMON)/kasmbl.o \
	$(COMMON)/sky.o 
OBJECTS2 = $(COMMON)/crossx.o
OBJECTS_BR = $(BRICK)/brshape.o $(BRICK)/brickb.o $(BRICK)/brcubic.o \
	$(BRICK)/brstress.o $(BRICK)/brwriter2.o 

PROGS = br2 

all: br2

br2: fembr2.o \
		brmemory.o \
		br2reader.o \
		br2id.o \
		brid.o \
		br2shape.o \
		brick2b.o \
		brcasmbl.o \
		br2conj.o \
		br2bound.o \
		brbound.o \
		br2kasmbl.o \
		br2writer.o \
		../brick/brconst.h \
		br2struct.h 
	$(CC) -o br2 fembr2.o \
		brmemory.o \
		br2reader.o \
		br2id.o \
		brid.o \
		br2shape.o \
		brick2b.o \
		brcasmbl.o \
		br2conj.o \
		br2bound.o \
		brbound.o \
		br2kasmbl.o \
		br2writer.o $(OBJECTS_BR) $(OBJECTS) $(LIBS)

br2writer.o: br2writer.c
	$(CC) -c br2writer.c 

br2kasmbl.o: br2kasmbl.c
	$(CC) -c br2kasmbl.c 

brick2b.o: brick2b.c
	$(CC) -c brick2b.c 

br2shape.o: br2shape.c
	$(CC) -c br2shape.c 

brbound.o: ../brick/brbound.c
	$(CC) -c ../brick/brbound.c -DBRICK2

br2bound.o: br2bound.c
	$(CC) -c br2bound.c 

br2conj.o: br2conj.c
	$(CC) -c br2conj.c 

brcasmbl.o: brcasmbl.c
	$(CC) -c brcasmbl.c 

brid.o: ../brick/brid.c
	$(CC) -c ../brick/brid.c -DBRICK2

br2id.o: br2id.c
	$(CC) -c br2id.c

br2reader.o: br2reader.c
	$(CC) -c br2reader.c

brmemory.o: ../brick/brmemory.c
	$(CC) -c ../brick/brmemory.c -DBRICK2

fembr2.o: fembr2.c
	$(CC) -c fembr2.c

clean:
	-rm *.o $(PROGS)
