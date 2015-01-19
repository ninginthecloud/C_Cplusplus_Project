#include "graph.h"

int main()
{
char graphfile[] = "data.txt";
char dotfile[] = "data.dot";
int nvertices = -1;
int** graph = readgraph(graphfile,nvertices);
//print out the graph in Graphviz format
printgraph(dotfile,graph,nvertices);
//find the connected component of vertex 1
findConComp(0,graph,nvertices);
//find the connected component of vertex 2
findConComp(1,graph,nvertices);
//find the connected component of vertex 5
findConComp(4,graph,nvertices);
//free memory
freegraph(graph,nvertices);
return(1);
}