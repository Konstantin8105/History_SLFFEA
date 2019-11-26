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

#May work on Sun
#C = /usr/lang/gcc  -g

#Dubugging Compile for Linux
CC = /usr/bin/gcc  -g

#Optimization Compile for Linux
#CC = /usr/bin/gcc  -O3

LIBS = -lm
COMMON = ../../common
BRICK = ../../brick/brick
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/dotx.o \
	$(COMMON)/lm.o $(COMMON)/idiag.o $(COMMON)/kasmbl.o  \
	$(COMMON)/sky.o $(COMMON)/qr.o $(COMMON)/eigen.o \
	$(COMMON)/eigenprint.o
OBJECTS_BR = $(BRICK)/brcubic.o

OBJECTS2 = $(COMMON)/crossx.o 

PROGS = shl shp shplane

all: brcubic.o shl shp shplane

shl: femsh.o \
		shmemory.o \
		shreader.o \
		shid.o \
		shshape.o \
		shellb.o \
		shkasmbl.o\
		shconj.o\
		shbound.o\
		shmasmbl.o \
		shlanczos.o \
		shmassp.o \
		shwriter.o\
		shconst.h \
		shstruct.h
	$(CC) -o shl femsh.o \
		shmemory.o \
		shreader.o \
		shid.o \
		shshape.o \
		shellb.o \
		shkasmbl.o\
		shconj.o\
		shbound.o\
		shmasmbl.o \
		shlanczos.o \
		shmassp.o \
		shwriter.o $(OBJECTS) $(OBJECTS_BR) $(OBJECTS2) $(LIBS)

brcubic.o:
	cd $(BRICK) ; make brcubic.o

shplane: shplane.o \
		shplane.o
	$(CC) -o shplane shplane.o $(LIBS) 

shp: shp.o \
		shpatch.o
	$(CC) -o shp shpatch.o $(LIBS) 

shplane.o: shplane.c
	$(CC) -c shplane.c

shp.o: shpatch.c
	$(CC) -c shpatch.c

shwriter.o: shwriter.c
	$(CC) -c shwriter.c 

shmassp.o: shmassp.c
	$(CC) -c shmassp.c 

shlanczos.o: shlanczos.c
	$(CC) -c shlanczos.c 

shmasmbl.o: shmasmbl.c
	$(CC) -c shmasmbl.c

#shstress.o: shstress.c
#	$(CC) -c shstress.c

shbound.o: shbound.c
	$(CC) -c shbound.c

shconj.o: shconj.c
	$(CC) -c shconj.c

shkasmbl.o: shkasmbl.c
	$(CC) -c shkasmbl.c

shellb.o: shellb.c
	$(CC) -c shellb.c 

shshape.o: shshape.c
	$(CC) -c shshape.c 

shid.o: shid.c
	$(CC) -c shid.c

shreader.o: shreader.c
	$(CC) -c shreader.c

shmemory.o: shmemory.c
	$(CC) -c shmemory.c

femsh.o: femsh.c
	$(CC) -c femsh.c

clean:
	-rm *.o $(PROGS)

