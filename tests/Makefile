
CC = g++
CFLAGS = -g 
INCL = -I../include/ -I/usr/local/bin/eigen
LIBS = --std=c++17

DEPS = ../lib/axis_convergence.o ../lib/ellipse_functions.o ../lib/helpers.o ../lib/constant_density_shells.o 
EXEC = test_basic test_cds test_ellipse_functions

all : $(EXEC)

%.o: %.cc $(DEPS) 
	$(CC) $(CFLAGS) $(INCL) $(LIBS) -c -o $@ $< 

test_basic : $(DEPS) test_basic.o
	$(CC) $(CFLAGS) $(INCL) $(LIBS) -o test_basic test_basic.o $(DEPS) 

test_ellipse_functions : $(DEPS) test_ellipse_functions.o
	$(CC) $(CFLAGS) $(INCL) $(LIBS) -o test_ellipse_functions test_ellipse_functions.o $(DEPS)

test_cds : $(DEPS) test_cds.o
	$(CC) $(CFLAGS) $(INCL) $(LIBS) -o test_cds test_cds.o $(DEPS)

clean:
	rm -f $(DEPS) $(EXEC) $(EXEC_O)


