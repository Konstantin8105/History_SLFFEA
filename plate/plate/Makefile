# SLFFEA Makefile
# Version:  1.3
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
	$(COMMON)/eigenprint.o

PROGS = pl plp

all: pl plp

pl: fempl.o \
		plmemory.o \
		plreader.o \
		plid.o \
		plshape.o \
		plateb.o \
		plkasmbl.o\
		plconj.o\
		plbound.o\
		plstress.o\
		plmasmbl.o \
		pllanczos.o \
		plmassp.o \
		plwriter.o\
		plarea.o\
		plconst.h \
		plstruct.h
	$(CC) -o pl fempl.o \
		plmemory.o \
		plreader.o \
		plid.o \
		plshape.o \
		plateb.o \
		plkasmbl.o\
		plconj.o\
		plbound.o\
		plstress.o\
		plmasmbl.o \
		pllanczos.o \
		plmassp.o \
		plwriter.o \
		plarea.o $(OBJECTS) $(LIBS)

plp: plp.o \
		plpatch.o
	$(CC) -o plp plpatch.o $(LIBS) 

plp.o: plpatch.c
	$(CC) -c plpatch.c

plarea.o: plarea.c
	$(CC) -c plarea.c 

plwriter.o: plwriter.c
	$(CC) -c plwriter.c 

plmassp.o: plmassp.c
	$(CC) -c plmassp.c

pllanczos.o: pllanczos.c
	$(CC) -c pllanczos.c

plmasmbl.o: plmasmbl.c
	$(CC) -c plmasmbl.c

plstress.o: plstress.c
	$(CC) -c plstress.c

plbound.o: plbound.c
	$(CC) -c plbound.c

plconj.o: plconj.c
	$(CC) -c plconj.c

plkasmbl.o: plkasmbl.c
	$(CC) -c plkasmbl.c

plateb.o: plateb.c
	$(CC) -c plateb.c 

plshape.o: plshape.c
	$(CC) -c plshape.c 

plid.o: plid.c
	$(CC) -c plid.c

plreader.o: plreader.c
	$(CC) -c plreader.c

plmemory.o: plmemory.c
	$(CC) -c plmemory.c

fempl.o: fempl.c
	$(CC) -c fempl.c

clean:
	-rm *.o $(PROGS)
