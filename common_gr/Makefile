# SLFFEA Makefile
# Version:  1.2
# Copyright (C) 1999, 2000, 2001  San Le (slffea.juno.com)
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


# CC = /usr/lang/acc -g
CC = /usr/bin/gcc

all: pwinston.o axes.o axeslabel.o mshcommon.o crtlcommon.o \
	printext.o screenshot.o filecheck.o scalcom.o \
	mytext.o pointlocate.o quat.o crossplane.o

pwinston.o: pwinston.c 
	$(CC) -c pwinston.c

axeslabel.o: axeslabel.c 
	$(CC) -c axeslabel.c

axes.o: axes.c 
	$(CC) -c axes.c

mshcommon.o: mshcommon.c 
	$(CC) -c mshcommon.c

crtlcommon.o: crtlcommon.c 
	$(CC) -c crtlcommon.c 

printext.o: printext.c 
	$(CC) -c printext.c

screenshot.o: screenshot.c
	$(CC) -c screenshot.c 

filecheck.o: filecheck.c
	$(CC) -c filecheck.c 

scalcom.o: scalcom.c
	$(CC) -c scalcom.c 

mytext.o: mytext.c
	$(CC) -c mytext.c 

pointlocate.o: pointlocate.c
	$(CC) -c pointlocate.c 

quat.o: quat.c
	$(CC) -c quat.c 

crossplane.o: crossplane.c
	$(CC) -c crossplane.c 

clean:
	-rm *.o 

