# SLFFEA Makefile
# Version:  1.1
# Copyright (C) 1999  San Le (slffea.juno.com)
# http://www.geocities.com/Athens/2099/slffea.html
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
	$(COMMON)/eigenprint.o $(COMMON)/dyadx.o

PROGS = qd3

all: qd3

qd3: femqd3.o \
		qd3memory.o \
		qd3reader.o \
		qd3id.o \
		qd3shape.o \
		quad3b.o \
		qd3kasmbl.o\
		qd3bound.o \
		qd3lanczos.o\
		qd3massp.o\
		qd3writer.o\
		qd3const.h \
		qd3struct.h
	$(CC) -o qd3 femqd3.o \
		qd3memory.o \
		qd3reader.o \
		qd3id.o \
		qd3shape.o \
		quad3b.o \
		qd3kasmbl.o\
		qd3bound.o \
		qd3lanczos.o\
		qd3massp.o\
		qd3writer.o $(OBJECTS) $(LIBS)


#qd3: femqd3.o \
#		qd3memory.o \
#		qd3reader.o \
#		qd3id.o \
#		qd3shape.o \
#		quad3b.o \
#		qd3kasmbl.o\
#		qd3conj.o\
#		qd3bound.o \
#		qd3stress.o\
#		qd3masmbl.o\
#		qd3lanczos.o\
#		qd3massp.o\
#		qd3writer.o\
#		qd3const.h \
#		qd3struct.h
#	$(CC) -o qd3 femqd3.o \
#		qd3memory.o \
#		qd3reader.o \
#		qd3id.o \
#		qd3shape.o \
#		quad3b.o \
#		qd3kasmbl.o\
#		qd3conj.o\
#		qd3bound.o \
#		qd3stress.o\
#		qd3masmbl.o\
#		qd3lanczos.o\
#		qd3massp.o\
#		qd3writer.o $(OBJECTS) $(LIBS)

qd3writer.o: qd3writer.c
	$(CC) -c qd3writer.c 

qd3massp.o: qd3massp.c
	$(CC) -c qd3massp.c 

qd3lanczos.o: qd3lanczos.c
	$(CC) -c qd3lanczos.c 

qd3masmbl.o: qd3masmbl.c
	$(CC) -c qd3masmbl.c 

qd3stress.o: qd3stress.c
	$(CC) -c qd3stress.c

qd3bound.o: qd3bound.c
	$(CC) -c qd3bound.c 

#qd3conj.o: qd3conj.c
#	$(CC) -c qd3conj.c

qd3kasmbl.o: qd3kasmbl.c
	$(CC) -c qd3kasmbl.c

quad3b.o: quad3b.c
	$(CC) -c quad3b.c 

qd3shape.o: qd3shape.c
	$(CC) -c qd3shape.c 

qd3id.o: qd3id.c
	$(CC) -c qd3id.c

qd3reader.o: qd3reader.c
	$(CC) -c qd3reader.c

qd3memory.o: qd3memory.c
	$(CC) -c qd3memory.c

femqd3.o: femqd3.c
	$(CC) -c femqd3.c

clean:
	-rm *.o $(PROGS)
