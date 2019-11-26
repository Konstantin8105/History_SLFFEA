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
TRI = ../../tri/tri
QUAD = ../../quad/quad
OBJECTS = $(COMMON)/matxt.o $(COMMON)/matx.o $(COMMON)/dotx.o \
	$(COMMON)/lm.o $(COMMON)/idiag.o $(COMMON)/kasmbl.o \
	$(COMMON)/sky.o $(COMMON)/qr.o $(COMMON)/eigen.o \
	$(COMMON)/eigenprint.o $(COMMON)/crossx.o \
	$(COMMON)/rotate4.o $(COMMON)/rotate2.o $(COMMON)/rotate.o \
	$(COMMON)/cubic.o $(COMMON)/local_vec.o $(COMMON)/local_vec2.o
OBJECTS_TRI = $(TRI)/trshape.o $(TRI)/trib.o 
OBJECTS_QD = $(QUAD)/quadb.o 

PROGS = pl plp

all: quadb.o trshape.o trib.o pl plp

pl: fempl.o \
		plmemory.o \
		plreader.o \
		plid.o \
		plshape.o \
		plateb.o \
		plkasmbl.o\
		plkasmbl2.o\
		plconj.o\
		plconj2.o\
		plbound.o\
		plstress.o\
		plmasmbl.o \
		plmasmbl2.o \
		pllanczos.o \
		plmassp.o \
		plmassp2.o \
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
		plkasmbl2.o\
		plconj.o\
		plconj2.o\
		plbound.o\
		plstress.o\
		plmasmbl.o \
		plmasmbl2.o \
		pllanczos.o \
		plmassp.o \
		plmassp2.o \
		plwriter.o \
		plarea.o $(OBJECTS) $(OBJECTS_QD) $(OBJECTS_TRI) $(LIBS)

trshape.o:
	cd $(TRI) ; make trshape.o 

trib.o:
	cd $(TRI) ; make trib.o 

quadb.o:
	cd $(QUAD) ; make quadb.o 

plp: plp.o \
		plpatch.o
	$(CC) -o plp plpatch.o $(LIBS) 

plp.o: plpatch.c
	$(CC) -c plpatch.c

plarea.o: plarea.c
	$(CC) -c plarea.c 

plwriter.o: plwriter.c
	$(CC) -c plwriter.c 

plmassp2.o: plmassp2.c
	$(CC) -c plmassp2.c

plmassp.o: plmassp.c
	$(CC) -c plmassp.c

pllanczos.o: pllanczos.c
	$(CC) -c pllanczos.c

plmasmbl2.o: plmasmbl2.c
	$(CC) -c plmasmbl2.c

plmasmbl.o: plmasmbl.c
	$(CC) -c plmasmbl.c

plstress.o: plstress.c
	$(CC) -c plstress.c

plbound.o: plbound.c
	$(CC) -c plbound.c

plconj2.o: plconj2.c
	$(CC) -c plconj2.c

plconj.o: plconj.c
	$(CC) -c plconj.c

plkasmbl2.o: plkasmbl2.c
	$(CC) -c plkasmbl2.c

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
