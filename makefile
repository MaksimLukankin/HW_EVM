CC=mpic++
CFLAGS= -O3 -fstack-protector-all -g  -Wall -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wno-suggest-attribute=format 

all: a.out

a.out: main.o io.o Jordan_mpi.o headers.h
	$(CC) $(CFLAGS) main.o io.o Jordan_mpi.o headers.h -o a.out
main.o: main.cpp
	$(CC) -c $(CFLAGS) main.cpp
io.o: io.cpp
	$(CC) -c $(CFLAGS) io.cpp
Jordan_mpi.o: Jordan_mpi.cpp
	$(CC) -c $(CFLAGS) Jordan_mpi.cpp
clean:
	rm -rf *.o a.out
