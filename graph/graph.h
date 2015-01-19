#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

int** allocgraph(int nvertices);
void freegraph(int**& graph,int nvertices);
int** readgraph(char* filename,int& nvertices);
void printgraph(char* filename,int** graph,int nvertices);
void findneighbor(int** graph,int* record, int nvertices, int myvertex);
void findConComp(int myvertex,int** graph,int nvertices);