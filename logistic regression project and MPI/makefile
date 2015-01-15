MPI_LSTDFLG = -lstdc++ -llapack -llapack_atlas -lcblas -lm -lstdc++ -lgsl -lgslcblas
MPI_INCLUDE = -I/usr/include/ -I/usr/include/atlas/
MPI_LIB = -L/usr/lib/atlas/ -L/usr/include/atlas_base/
MPI_OBJS = parallelBayes

all:	${MPI_OBJS}
	rm -f *.o

matrices.o: matrices.cpp matrices.h
	mpic++ -g -c matrices.cpp -o matrices.o ${MPI_INCLUDE}

Bayes.o: Bayes.cpp Bayes.h
	mpic++ -g -c Bayes.cpp -o Bayes.o ${MPI_INCLUDE}

main.o: main.cpp matrices.h
	mpic++ -g -c main.cpp -o main.o ${MPI_INCLUDE}

parallelBayes: main.o matrices.o Bayes.o
	mpic++ main.o Bayes.o matrices.o -o parallelBayes ${MPI_LIB} ${MPI_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MPI_OBJS}