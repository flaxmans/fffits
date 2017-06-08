#
TARGET=fffits
CC=gcc
CFLAGS=-lm -I/usr/local/include/gsl -L/usr/local/lib -lgsl -lgslcblas -O3

# rule to build fffits
$(TARGET) : $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).c

# rule to clean
clean :
	rm -f $(TARGET)


