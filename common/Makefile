# SLFFEA Makefile
# Version:  1.5
# Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005  San Le (slffea@yahoo.com)
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
	qr.o eigen.o eigenprint.o rotate.o rotate2.o rotate3.o \
	rotate4.o dyadx.o cubic.o nodenorm.o local_vec.o local_vec2.o \
	gnuplotwrite2.o compare

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

rotate2.o: rotate2.c
	$(CC) -c rotate2.c 

rotate3.o: rotate3.c
	$(CC) -c rotate3.c 

rotate4.o: rotate4.c
	$(CC) -c rotate4.c 

cubic.o: cubic.c
	$(CC) -c cubic.c 

nodenorm.o: nodenorm.c
	$(CC) -c nodenorm.c 

local_vec.o: local_vec.c
	$(CC) -c local_vec.c 

local_vec2.o: local_vec2.c
	$(CC) -c local_vec2.c 

gnuplotwrite2.o: gnuplotwrite2.c
	$(CC) -c gnuplotwrite2.c 

compare: compare.o \
		compare.o
	$(CC) -o compare compare.o $(LIBS)

compare.o: compare.c
	$(CC) -c compare.c 

clean:
	-rm *.o 
