GRAPHS_LSTDFLG = -lstdc++ -llapack -lcblas -lm
GRAPHS_INCLUDE = -I/usr/include/
GRAPHS_LIB = -L/usr/lib/atlas/
GRAPHS_OBJS = graphs

all:	${GRAPHS_OBJS}
	rm -f *.o

graph.o: graph.cpp graph.h
	gcc -g -c graph.cpp -o graph.o ${GRAPHS_INCLUDE}

main.o: main.cpp graph.h
	gcc -g -c main.cpp -o main.o

graphs: main.o graph.o
	gcc main.o graph.o -o graphs ${GRAPHS_LIB} ${GRAPHS_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${GRAPHS_OBJS}