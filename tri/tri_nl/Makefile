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
TRI = ../tri
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/crossx.o \
	$(COMMON)/dotx.o $(COMMON)/cubic.o $(COMMON)/kasmbl.o \
	$(COMMON)/rotate2.o $(COMMON)/rotate3.o $(COMMON)/local_vec2.o \
	$(COMMON)/gnuplotwrite2.o
OBJECTS2 = $(COMMON)/crossx.o 
OBJECTS_TR = $(TRI)/memory2.o $(TRI)/trreader.o $(TRI)/id.o $(TRI)/trshape.o \
	$(TRI)/trib.o $(TRI)/bound.o $(TRI)/trwriter.o \
	$(TRI)/trconj.o $(TRI)/trkasmbl.o

#PROGS = ntr pvv
PROGS = ntr

#all: ntr pvv
all: ntr

ntr: femntr.o \
		trfmasmbl.o \
		trpasmbl2.o \
		../tri/trconst.h \
		../tri/trstruct.h 
	$(CC) -o ntr femntr.o \
		trfmasmbl.o \
		trpasmbl2.o $(OBJECTS_TR) $(OBJECTS) $(LIBS)

trpasmbl2.o: trpasmbl2.c
	$(CC) -c trpasmbl2.c 

trfmasmbl.o: trfmasmbl.c
	$(CC) -c trfmasmbl.c 

femntr.o: femntr.c
	$(CC) -c femntr.c

clean:
	-rm *.o $(PROGS)
