# SLFFEA Makefile
# Version:  1.4
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
BEAM = ../../beam/beam
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/dotx.o \
	$(COMMON)/lm.o $(COMMON)/idiag.o $(COMMON)/kasmbl.o \
	$(COMMON)/sky.o $(COMMON)/qr.o $(COMMON)/eigen.o \
	$(COMMON)/eigenprint.o $(COMMON)/rotate.o
PROGS = ts 

all: ts 

ts: femtruss.o \
		tsmemory.o \
		tsreader.o \
		tslocal_vec.o \
		id.o \
		tskasmbl.o \
		tsconj.o \
		bound.o \
		tsmasmbl.o \
		tslanczos.o \
		tsmassp.o \
		tswriter.o \
		tslength.o \
		tsconst.h \
		tsstruct.h 
	$(CC) -o ts femtruss.o \
		tsmemory.o \
		tsreader.o \
		tslocal_vec.o \
		id.o \
		tskasmbl.o \
		tsconj.o \
		bound.o \
		tsmasmbl.o \
		tslanczos.o \
		tsmassp.o \
		tswriter.o \
		tslength.o $(OBJECTS) $(LIBS)


tslength.o: tslength.c
	$(CC) -c tslength.c 

tswriter.o: tswriter.c
	$(CC) -c tswriter.c 

tsmassp.o: tsmassp.c
	$(CC) -c tsmassp.c

tslanczos.o: tslanczos.c
	$(CC) -c tslanczos.c

tsmasmbl.o: tsmasmbl.c
	$(CC) -c tsmasmbl.c

bound.o: $(COMMON)/bound.c
	$(CC) -c $(COMMON)/bound.c -DTRUSS1

tsconj.o: tsconj.c
	$(CC) -c tsconj.c 

tskasmbl.o: tskasmbl.c
	$(CC) -c tskasmbl.c 

id.o: $(COMMON)/id.c
	$(CC) -c $(COMMON)/id.c -DTRUSS1

tslocal_vec.o: tslocal_vec.c
	$(CC) -c tslocal_vec.c

tsreader.o: tsreader.c
	$(CC) -c tsreader.c

tsmemory.o: tsmemory.c
	$(CC) -c tsmemory.c

femtruss.o: femtruss.c
	$(CC) -c femtruss.c

clean:
	-rm *.o $(PROGS)
