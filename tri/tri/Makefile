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
	$(COMMON)/sky.o $(COMMON)/qr.o $(COMMON)/eigen.o \
	$(COMMON)/eigenprint.o $(COMMON)/crossx.o \
	$(COMMON)/rotate2.o $(COMMON)/local_vec2.o

PROGS = tri trp

all: tri trp

tri: femtr.o \
		memory2.o \
		trreader.o \
		id.o \
		trshape.o \
		trib.o \
		trkasmbl.o\
		trconj.o\
		bound.o \
		trmasmbl.o\
		trlanczos.o\
		trmassp.o\
		trwriter.o\
		trconst.h \
		trstruct.h
	$(CC) -o tri femtr.o \
		memory2.o \
		trreader.o \
		id.o \
		trshape.o \
		trib.o \
		trkasmbl.o\
		trconj.o\
		bound.o \
		trmasmbl.o\
		trlanczos.o\
		trmassp.o\
		trwriter.o $(OBJECTS) $(LIBS)

trp: trp.o \
		trpatch.o
	$(CC) -o trp trpatch.o $(LIBS) 

trp.o: trpatch.c
	$(CC) -c trpatch.c

trwriter.o: trwriter.c
	$(CC) -c trwriter.c 

trmassp.o: trmassp.c
	$(CC) -c trmassp.c 

trlanczos.o: trlanczos.c
	$(CC) -c trlanczos.c 

trmasmbl.o: trmasmbl.c
	$(CC) -c trmasmbl.c 

bound.o: $(COMMON)/bound.c
	$(CC) -c $(COMMON)/bound.c -DTRI1

trconj.o: trconj.c
	$(CC) -c trconj.c

trkasmbl.o: trkasmbl.c
	$(CC) -c trkasmbl.c

trib.o: trib.c
	$(CC) -c trib.c 

trshape.o: trshape.c
	$(CC) -c trshape.c 

id.o: $(COMMON)/id.c
	$(CC) -c $(COMMON)/id.c -DTRI1

trreader.o: trreader.c
	$(CC) -c trreader.c

memory2.o: $(COMMON)/memory2.c
	$(CC) -c $(COMMON)/memory2.c -DTRI1

femtr.o: femtr.c
	$(CC) -c femtr.c

clean:
	-rm *.o $(PROGS)
