
CC = g++
CFLAGS = -g -Wall
INCL = -I../include/
LIBS = --std=c++17

DEPS = ../lib/ellipse_functions.o ../lib/helpers.o
EXEC = test_ellipse_functions

all : $(EXEC)

%.o: %.cc $(DEPS) 
	$(CC) $(CFLAGS) $(INCL) $(LIBS) -c -o $@ $< 

test_ellipse_functions : $(DEPS) test_ellipse_functions.o
	$(CC) $(CFLAGS) $(INCL) $(LIBS) -o test_ellipse_functions test_ellipse_functions.o $(DEPS)

clean:
	rm -f $(DEPS) $(EXEC) test_ellipse_functions.o 


