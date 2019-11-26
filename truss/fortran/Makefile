F77 = /usr/bin/g77  -g
LIBS = -lm
PROGS = ftf 

all: ftf 

ftf: truss.o \
		tsconst.hf 
	$(F77) -o ftf truss.o $(LIBS) 

truss.o: truss.f
	$(F77) -c truss.f 

clean:
	-rm *.o $(PROGS)
