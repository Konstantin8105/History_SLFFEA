# SLFFEA Makefile
# Version:  1.3
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
PROGS = we wep

all: we wep

we: femwe.o \
		wememory.o \
		wereader.o \
		weid.o \
		weshape.o \
		wedgeb.o \
		wekasmbl.o \
		weconj.o \
		webound.o \
		westress.o \
		wemasmbl.o \
		welanczos.o \
		wemassp.o \
		wewriter.o \
		wevolume.o \
		wewriter2.o \
		weconst.h \
		westruct.h 
	$(CC) -o we femwe.o \
		wememory.o \
		wereader.o \
		weid.o \
		weshape.o \
		wedgeb.o \
		wekasmbl.o \
		weconj.o \
		webound.o \
		westress.o \
		wemasmbl.o \
		welanczos.o \
		wemassp.o \
		wewriter.o \
		wevolume.o \
		wewriter2.o $(OBJECTS) $(LIBS)

wep: wep.o \
		wepatch.o
	$(CC) -o wep wepatch.o $(OBJECTS2) $(LIBS)
		
wep.o: wepatch.c
	$(CC) -c wepatch.c

wewriter2.o: wewriter2.c
	$(CC) -c wewriter2.c 

wevolume.o: wevolume.c
	$(CC) -c wevolume.c 

wewriter.o: wewriter.c
	$(CC) -c wewriter.c 

wemassp.o: wemassp.c
	$(CC) -c wemassp.c 

welanczos.o: welanczos.c
	$(CC) -c welanczos.c 

wemasmbl.o: wemasmbl.c
	$(CC) -c wemasmbl.c 

westress.o: westress.c
	$(CC) -c westress.c 

webound.o: webound.c
	$(CC) -c webound.c -DWEDGE1

weconj.o: weconj.c
	$(CC) -c weconj.c 

wekasmbl.o: wekasmbl.c
	$(CC) -c wekasmbl.c 

wedgeb.o: wedgeb.c
	$(CC) -c wedgeb.c 

weshape.o: weshape.c
	$(CC) -c weshape.c 

weid.o: weid.c
	$(CC) -c weid.c -DWEDGE1

wereader.o: wereader.c
	$(CC) -c wereader.c

wememory.o: wememory.c
	$(CC) -c wememory.c -DWEDGE1

femwe.o: femwe.c
	$(CC) -c femwe.c

clean:
	-rm *.o $(PROGS)
