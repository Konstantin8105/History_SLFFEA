# SLFFEA Makefile
# Version:  1.3
# Copyright (C) 1999, 2000  San Le (slffea@yahoo.com)
# http://www.geocities.com/slffea/slffea.html
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
TRUSS = ../truss
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/dotx.o \
	$(COMMON)/rotate.o
OBJECTS_TS = $(TRUSS)/tsmemory.o $(TRUSS)/tsreader.o $(TRUSS)/id.o \
	$(TRUSS)/bound.o $(TRUSS)/tswriter.o $(TRUSS)/tslength.o

PROGS = nts 

all: nts 

nts: femntruss.o \
		tsfmasmbl.o \
		tspasmbl2.o \
		../truss/tsconst.h \
		../truss/tsstruct.h 
	$(CC) -o nts femntruss.o \
		tsfmasmbl.o \
		tspasmbl2.o $(OBJECTS_TS) $(OBJECTS) $(LIBS)


tspasmbl2.o: tspasmbl2.c
	$(CC) -c tspasmbl2.c

tsfmasmbl.o: tsfmasmbl.c
	$(CC) -c tsfmasmbl.c

femntruss.o: femntruss.c
	$(CC) -c femntruss.c

clean:
	-rm *.o $(PROGS)
