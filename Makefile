#
TARGET=fffits
CC=gcc
CFLAGS=-g -DDEBUG
GSLFLAGS=-I/usr/local/include/gsl
LIBFLAGS=-lm -lgsl -lgslcblas -L/usr/local/lib

# rule to build fffits
$(TARGET) : $(TARGET).o initialization.o dataRecording.o
	$(CC) $(CFLAGS) $(GSLFLAGS) $(LIBFLAGS) -o $(TARGET) initialization.o dataRecording.o $(TARGET).o

# building a .o from a .c
.c.o:
	$(CC) $(CFLAGS) $(GSLFLAGS) -c $<

# rule to clean
clean :
	rm -f $(TARGET)
	rm -f $(TARGET).o
	rm -f initialization.o
	rm -f dataRecording.o

