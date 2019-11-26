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
BRICK = ../brick
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/crossx.o \
	$(COMMON)/dotx.o $(COMMON)/cubic.o 
OBJECTS2 = $(COMMON)/crossx.o 
OBJECTS_BR = $(BRICK)/memory.o $(BRICK)/brreader.o $(BRICK)/id.o $(BRICK)/brshape.o \
	$(BRICK)/brickb.o $(BRICK)/bound.o $(BRICK)/brwriter.o $(BRICK)/brvolume.o  \
	$(BRICK)/brwriter2.o $(BRICK)/brconj.o 

#PROGS = nbr pvv
PROGS = nbr

#all: nbr pvv
all: nbr

nbr: femnbr.o \
		brfmasmbl.o \
		brpasmbl.o \
		brpasmbl2.o \
		../brick/brconst.h \
		../brick/brstruct.h 
	$(CC) -o nbr femnbr.o \
		brfmasmbl.o \
		brpasmbl.o \
		brpasmbl2.o $(OBJECTS_BR) $(OBJECTS) $(LIBS)

brpasmbl2.o: brpasmbl2.c
	$(CC) -c brpasmbl2.c 

brpasmbl.o: brpasmbl.c
	$(CC) -c brpasmbl.c 

brfmasmbl.o: brfmasmbl.c
	$(CC) -c brfmasmbl.c 

femnbr.o: femnbr.c
	$(CC) -c femnbr.c

clean:
	-rm *.o $(PROGS)
