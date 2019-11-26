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


all: wepost

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_WE = ../wedge/wememory.o ../wedge/wereader.o 

OBJECTS_WEH = ../wedge/weconst.h ../wedge/westruct.h 

COMMON = ../../common

COMMON_GR = ../../common_gr

OBJECTS_GR = $(COMMON_GR)/pwinston.o $(COMMON_GR)/axes.o \
	$(COMMON_GR)/mshcommon.o $(COMMON_GR)/crtlcommon.o \
	$(COMMON_GR)/printext.o $(COMMON_GR)/screenshot.o \
	$(COMMON_GR)/filecheck.o $(COMMON_GR)/pointlocate.o \
	$(COMMON_GR)/axeslabel.o $(COMMON_GR)/quat.o \
	$(COMMON_GR)/crossplane.o

OBJECTS_GRH = $(COMMON_GR)/color_gr.h $(COMMON_GR)/control.h 

OBJECTS2 = $(COMMON)/crossx.o $(COMMON)/matx.o

PROGS = wepost

CC = /usr/bin/gcc -g
LIBS = -lm

wepost: wepost.o \
		wememgr.o \
		wereadgr.o \
		wereadgr2.o \
		weparam.o \
		weparam2.o \
		weset.o \
		wenorm_vec.o \
		webc_vec.o \
		wemshdsp.o \
		wemesh.o \
		wemshkey.o \
		wemenu.o \
		wegetnew.o \
		wecrtldsp.o \
		wecrtlmse.o \
		werender.o \
		$(OBJECTS_WEH) \
		$(OBJECTS_GRH) \
		wegui.h \
		westrcgr.h 
	$(CC) -o wepost $(LIB1) wepost.o \
		wememgr.o \
		wereadgr.o \
		wereadgr2.o \
		weparam.o \
		weparam2.o \
		weset.o \
		wenorm_vec.o \
		webc_vec.o \
		wemshdsp.o \
		wemesh.o \
		wemshkey.o \
		wemenu.o \
		wegetnew.o \
		wecrtldsp.o \
		wecrtlmse.o \
		werender.o \
		$(OBJECTS_WE) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

werender.o: werender.c 
	$(CC) -c werender.c -DWEDGE1

wecrtlmse.o: wecrtlmse.c 
	$(CC) -c wecrtlmse.c 

wecrtldsp.o: wecrtldsp.c 
	$(CC) -c wecrtldsp.c -DWEDGE1

wegetnew.o: wegetnew.c 
	$(CC) -c wegetnew.c 

wemenu.o: wemenu.c 
	$(CC) -c wemenu.c 

wemshkey.o: wemshkey.c 
	$(CC) -c wemshkey.c 

wemesh.o: wemesh.c 
	$(CC) -c wemesh.c -DWEDGE1

wemshdsp.o: wemshdsp.c 
	$(CC) -c wemshdsp.c -DWEDGE1

webc_vec.o: webc_vec.c 
	$(CC) -c webc_vec.c -DWEDGE1

wenorm_vec.o: wenorm_vec.c 
	$(CC) -c wenorm_vec.c 

weset.o: weset.c 
	$(CC) -c weset.c 

weparam2.o: weparam2.c 
	$(CC) -c weparam2.c 

weparam.o: weparam.c 
	$(CC) -c weparam.c 

wereadgr2.o: wereadgr2.c 
	$(CC) -c wereadgr2.c -DWEDGE1

wereadgr.o: wereadgr.c 
	$(CC) -c wereadgr.c -DWEDGE1

wememgr.o: wememgr.c 
	$(CC) -c wememgr.c -DWEDGE1

wepost.o: wepost.c 
	$(CC) -c wepost.c

clean:
	-rm $(PROGS)
	-rm *.o 
