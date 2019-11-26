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


# CC = /usr/lang/acc -g
CC = /usr/bin/gcc

all: pwinston.o axes.o mshcommon.o crtlcommon.o \
	printext.o screenshot.o filecheck.o

pwinston.o: pwinston.c 
	$(CC) -c pwinston.c

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

clean:
	-rm *.o 

