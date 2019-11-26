# SLFFEA Makefile
# Version:  1.4
# Copyright (C) 1999, 2000, 2001, 2002  San Le (slffea@yahoo.com)
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

PROGS = te tep

all: te tep

te: femte.o \
		memory2.o \
		tereader.o \
		id.o \
		teshape.o \
		tetrab.o \
		tekasmbl.o \
		teconj.o \
		bound.o \
		temasmbl.o \
		telanczos.o \
		temassp.o \
		tewriter.o \
		tewriter2.o \
		teconst.h \
		testruct.h 
	$(CC) -o te femte.o \
		memory2.o \
		tereader.o \
		id.o \
		teshape.o \
		tetrab.o \
		tekasmbl.o \
		teconj.o \
		bound.o \
		temasmbl.o \
		telanczos.o \
		temassp.o \
		tewriter.o \
		tewriter2.o $(OBJECTS) $(LIBS)

tep: tep.o \
		tepatch.o
	$(CC) -o tep tepatch.o $(OBJECTS2) $(LIBS)
		
tep.o: tepatch.c
	$(CC) -c tepatch.c

tewriter2.o: tewriter2.c
	$(CC) -c tewriter2.c 

tewriter.o: tewriter.c
	$(CC) -c tewriter.c 

temassp.o: temassp.c
	$(CC) -c temassp.c 

telanczos.o: telanczos.c
	$(CC) -c telanczos.c 

temasmbl.o: temasmbl.c
	$(CC) -c temasmbl.c 

bound.o: $(COMMON)/bound.c
	$(CC) -c $(COMMON)/bound.c -DTETRA1

teconj.o: teconj.c
	$(CC) -c teconj.c 

tekasmbl.o: tekasmbl.c
	$(CC) -c tekasmbl.c 

tetrab.o: tetrab.c
	$(CC) -c tetrab.c 

teshape.o: teshape.c
	$(CC) -c teshape.c 

id.o: $(COMMON)/id.c
	$(CC) -c $(COMMON)/id.c -DTETRA1

tereader.o: tereader.c
	$(CC) -c tereader.c

memory2.o: $(COMMON)/memory2.c
	$(CC) -c $(COMMON)/memory2.c -DTETRA1

femte.o: femte.c
	$(CC) -c femte.c

clean:
	-rm *.o $(PROGS)
