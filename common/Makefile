# SLFFEA Makefile
# Version:  1.1
# Copyright (C) 1999  San Le (slffea.juno.com)
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

all: crossx.o dotx.o matx.o matxt.o lm.o idiag.o kasmbl.o sky.o \
	qr.o eigen.o eigenprint.o rotate.o dyadx.o

crossx.o: crossx.c
	$(CC) -c crossx.c 

dotx.o: dotx.c
	$(CC) -c dotx.c 

dyadx.o: dyadx.c
	$(CC) -c dyadx.c 

matx.o: matx.c
	$(CC) -c matx.c 

matxt.o: matxt.c
	$(CC) -c matxt.c 

idiag.o: idiag.c
	$(CC) -c idiag.c

lm.o: lm.c
	$(CC) -c lm.c

kasmbl.o: kasmbl.c
	$(CC) -c kasmbl.c 

sky.o: sky.c
	$(CC) -c sky.c 

qr.o: qr.c
	$(CC) -c qr.c 

eigen.o: eigen.c
	$(CC) -c eigen.c 

eigenprint.o: eigenprint.c
	$(CC) -c eigenprint.c 

rotate.o: rotate.c
	$(CC) -c rotate.c 

clean:
	-rm *.o 
