# SLFFEA Makefile
# Version:  1.2
# Copyright (C) 2008  San Le (slffea@yahoo.com)
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
TETRA = ../tetra
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/crossx.o \
	$(COMMON)/dotx.o $(COMMON)/cubic.o $(COMMON)/kasmbl.o \
	$(COMMON)/gnuplotwrite2.o
OBJECTS2 = $(COMMON)/crossx.o 
OBJECTS_TE = $(TETRA)/memory2.o $(TETRA)/tereader.o $(TETRA)/id.o $(TETRA)/teshape.o \
	$(TETRA)/tetrab.o $(TETRA)/bound.o $(TETRA)/tewriter.o \
	$(TETRA)/tewriter2.o $(TETRA)/teconj.o $(TETRA)/tekasmbl.o

#PROGS = nte pvv
PROGS = nte

#all: nte pvv
all: nte

nte: femnte.o \
		tefmasmbl.o \
		tepasmbl2.o \
		../tetra/teconst.h \
		../tetra/testruct.h 
	$(CC) -o nte femnte.o \
		tefmasmbl.o \
		tepasmbl2.o $(OBJECTS_TE) $(OBJECTS) $(LIBS)

tepasmbl2.o: tepasmbl2.c
	$(CC) -c tepasmbl2.c 

tefmasmbl.o: tefmasmbl.c
	$(CC) -c tefmasmbl.c 

femnte.o: femnte.c
	$(CC) -c femnte.c

clean:
	-rm *.o $(PROGS)
