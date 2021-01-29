# Makefile for compiling fffits
# To compile normally, just type "make" on the command line (but not the quotes)
# To compile in debugging mode, type "make DEBUG_MODE=Y" (but not the quotes)

TARGET=fffits
CC=gcc
GSLFLAGS=-I/usr/local/include/gsl
LIBFLAGS=-lm -lgsl -lgslcblas -L/usr/local/lib
PCGLIB=pcg-c-master/src/libpcg_random.a

# conditional compiling:
DEBUG_MODE?=no
ifeq "$(DEBUG_MODE)" "Y"
	CFLAGS=-g -DDEBUG
else 
	CFLAGS=-O3
endif

# rule to build fffits
$(TARGET) : $(TARGET).o initialization.o dataRecording.o pcg-c-master/src/libpcg_random.a
	$(CC) $(CFLAGS) $(GSLFLAGS) $(LIBFLAGS) -o $(TARGET) initialization.o dataRecording.o $(PCGLIB) $(TARGET).o

# building a .o from a .c
.c.o:
	$(CC) $(CFLAGS) $(GSLFLAGS) -c $<

# building pcg:
$(PCGLIB):
	cd pcg-c-master; make; cd ..

# rule to clean
clean :
	rm -f $(TARGET)
	rm -f $(TARGET).o
	rm -f initialization.o
	rm -f dataRecording.o
	rm -f dataRecording.o

cleanpcg : 
	cd pcg-c-master; make clean; cd ..

