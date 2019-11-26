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


all: plpost

LIB1 = -I/usr/local/include -O2 -funroll-loops -ansi -pedantic \
	 -ffast-math -I/usr/X11R6/include -DSHM 

LIB2 = -L/usr/local/lib -lglut -lGLU -lGL -lm -L/usr/X11/lib \
	-L/usr/X11R6/lib -lX11 -lXext -lXmu -lXt -lXi -lSM -lICE

OBJECTS_PL = ../plate/plmemory.o ../plate/plreader.o 

OBJECTS_PLH = ../plate/plconst.h ../plate/plstruct.h 

COMMON = ../../common

COMMON_GR = ../../common_gr

QUAD_GR = ../../quad/quad_gr

TRI_GR = ../../tri/tri_gr

OBJECTS_GR = $(COMMON_GR)/pwinston.o $(COMMON_GR)/axes.o \
	$(COMMON_GR)/mshcommon.o $(COMMON_GR)/crtlcommon.o \
	$(COMMON_GR)/printext.o $(COMMON_GR)/screenshot.o \
	$(COMMON_GR)/filecheck.o $(COMMON_GR)/pointlocate.o \
	$(COMMON_GR)/axeslabel.o $(COMMON_GR)/quat.o \
	$(COMMON_GR)/crossplane.o


OBJECTS_GRH = $(COMMON_GR)/color_gr.h $(COMMON_GR)/control.h 

OBJECTS2 = $(COMMON)/crossx.o $(COMMON)/matx.o

COMMON_GR = ../../common_gr

PROGS = plpost

CC = /usr/bin/gcc -g

plpost: plpost.o \
		plmemgr.o \
		plreadgr.o \
		plparam.o \
		plregetparam.o \
		plset.o \
		qdnorm_vec.o \
		trnorm_vec.o \
		plbc_vec.o \
		plmshdsp.o \
		plmesh.o \
		plmshkey.o \
		plmshkey2.o \
		plmenu.o \
		plgetnew.o \
		plcrtldsp.o \
		plcrtlmse.o \
		plrender.o \
		$(OBJECTS_PLH) \
		$(OBJECTS_GRH) \
		plgui.h \
		plstrcgr.h 
	$(CC) -o plpost $(LIB1) plpost.o \
		plmemgr.o \
		plreadgr.o \
		plparam.o \
		plregetparam.o \
		plset.o \
		qdnorm_vec.o \
		trnorm_vec.o \
		plbc_vec.o \
		plmshdsp.o \
		plmesh.o \
		plmshkey.o \
		plmshkey2.o \
		plmenu.o \
		plgetnew.o \
		plcrtldsp.o \
		plcrtlmse.o \
		plrender.o \
		$(OBJECTS_PL) $(OBJECTS_GR) $(OBJECTS2) $(LIB2) 

plrender.o: plrender.c 
	$(CC) -c plrender.c 

plcrtlmse.o: plcrtlmse.c 
	$(CC) -c plcrtlmse.c 

plcrtldsp.o: plcrtldsp.c 
	$(CC) -c plcrtldsp.c 

plgetnew.o: plgetnew.c 
	$(CC) -c plgetnew.c 

plmenu.o: plmenu.c 
	$(CC) -c plmenu.c 

plmshkey2.o: plmshkey2.c 
	$(CC) -c plmshkey2.c 

plmshkey.o: plmshkey.c 
	$(CC) -c plmshkey.c 

plmesh.o: plmesh.c 
	$(CC) -c plmesh.c 

plmshdsp.o: plmshdsp.c 
	$(CC) -c plmshdsp.c 

plbc_vec.o: plbc_vec.c 
	$(CC) -c plbc_vec.c 

qdnorm_vec.o: $(QUAD_GR)/qdnorm_vec.c 
	$(CC) -c $(QUAD_GR)/qdnorm_vec.c -DPLATE1

trnorm_vec.o: $(TRI_GR)/trnorm_vec.c 
	$(CC) -c $(TRI_GR)/trnorm_vec.c -DPLATE1

plset.o: plset.c 
	$(CC) -c plset.c 

plregetparam.o: plregetparam.c 
	$(CC) -c plregetparam.c 

plparam.o: plparam.c 
	$(CC) -c plparam.c 

plreadgr.o: plreadgr.c 
	$(CC) -c plreadgr.c 

plmemgr.o: plmemgr.c 
	$(CC) -c plmemgr.c 

plpost.o: plpost.c 
	$(CC) -c plpost.c

clean:
	-rm $(PROGS)
	-rm *.o 

