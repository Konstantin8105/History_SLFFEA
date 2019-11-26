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
	$(COMMON)/sky.o $(COMMON)/qr.o $(COMMON)/eigen.o \
	$(COMMON)/eigenprint.o $(COMMON)/cubic.o
OBJECTS2 = $(COMMON)/crossx.o
PROGS = br brp

all: br brp

br: fembr.o \
		brmemory.o \
		brreader.o \
		brid.o \
		brshape.o \
		brickb.o \
		brkasmbl.o \
		brconj.o \
		brbound.o \
		brstress.o \
		brmasmbl.o \
		brlanczos.o \
		brmassp.o \
		brwriter.o \
		brvolume.o \
		brwriter2.o \
		brconst.h \
		brstruct.h 
	$(CC) -o br fembr.o \
		brmemory.o \
		brreader.o \
		brid.o \
		brshape.o \
		brickb.o \
		brkasmbl.o \
		brconj.o \
		brbound.o \
		brstress.o \
		brmasmbl.o \
		brlanczos.o \
		brmassp.o \
		brwriter.o \
		brvolume.o \
		brwriter2.o $(OBJECTS) $(LIBS)

brp: brp.o \
		brpatch.o
	$(CC) -o brp brpatch.o $(OBJECTS2) $(LIBS)
		
brp.o: brpatch.c
	$(CC) -c brpatch.c

brwriter2.o: brwriter2.c
	$(CC) -c brwriter2.c 

brvolume.o: brvolume.c
	$(CC) -c brvolume.c 

brwriter.o: brwriter.c
	$(CC) -c brwriter.c 

brmassp.o: brmassp.c
	$(CC) -c brmassp.c 

brlanczos.o: brlanczos.c
	$(CC) -c brlanczos.c 

brmasmbl.o: brmasmbl.c
	$(CC) -c brmasmbl.c 

brstress.o: brstress.c
	$(CC) -c brstress.c 

brbound.o: brbound.c
	$(CC) -c brbound.c -DBRICK1

brconj.o: brconj.c
	$(CC) -c brconj.c 

brkasmbl.o: brkasmbl.c
	$(CC) -c brkasmbl.c 

brickb.o: brickb.c
	$(CC) -c brickb.c 

brshape.o: brshape.c
	$(CC) -c brshape.c 

brid.o: brid.c
	$(CC) -c brid.c -DBRICK1

brreader.o: brreader.c
	$(CC) -c brreader.c

brmemory.o: brmemory.c
	$(CC) -c brmemory.c -DBRICK1

fembr.o: fembr.c
	$(CC) -c fembr.c

clean:
	-rm *.o $(PROGS)
