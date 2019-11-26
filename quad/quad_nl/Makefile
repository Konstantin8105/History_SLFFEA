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
QUAD = ../quad
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/crossx.o \
	$(COMMON)/dotx.o $(COMMON)/cubic.o $(COMMON)/kasmbl.o \
	$(COMMON)/rotate2.o $(COMMON)/rotate3.o $(COMMON)/local_vec.o \
	$(COMMON)/gnuplotwrite2.o
OBJECTS2 = $(COMMON)/crossx.o 
OBJECTS_QD = $(QUAD)/memory.o $(QUAD)/qdreader.o $(QUAD)/id.o $(QUAD)/qdshape.o \
	$(QUAD)/quadb.o $(QUAD)/bound.o $(QUAD)/qdwriter.o $(QUAD)/qdarea.o  \
	$(QUAD)/qdconj.o $(QUAD)/qdkasmbl.o $(QUAD)/qdstress.o

#PROGS = nqd pvv
PROGS = nqd

#all: nqd pvv
all: nqd

nqd: femnqd.o \
		qdfmasmbl.o \
		qdpasmbl2.o \
		../quad/qdconst.h \
		../quad/qdstruct.h 
	$(CC) -o nqd femnqd.o \
		qdfmasmbl.o \
		qdpasmbl2.o $(OBJECTS_QD) $(OBJECTS) $(LIBS)

qdpasmbl2.o: qdpasmbl2.c
	$(CC) -c qdpasmbl2.c 

qdfmasmbl.o: qdfmasmbl.c
	$(CC) -c qdfmasmbl.c 

femnqd.o: femnqd.c
	$(CC) -c femnqd.c

clean:
	-rm *.o $(PROGS)
