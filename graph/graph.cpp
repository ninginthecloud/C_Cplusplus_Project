#include "graph.h"

//the graph is stored as an incidence matrix
//An incidence matrix is a symmetric square matrix
//with dimension equal to the number of vertices in the graph
//The (i,j) entry in this matrix takes value one
//if there is an edge between vertices i and j,
// and takes value zero otherwise.
int** allocgraph(int nvertices)
{
   int i;
   int** graph = new int*[nvertices];

   for(i=0;i<nvertices;i++)
   {
      graph[i] = new int[nvertices];
      //initially the graph has no edges
      memset(graph[i],0,nvertices*sizeof(int));
   }

   return(graph);
}

//free the memory for an incidence matrix
void freegraph(int**& graph,int nvertices)
{
   int i;

   for(i=0;i<nvertices;i++)
   {
      delete[] graph[i];
      graph[i] = NULL;
   }
   delete[] graph;
   graph = NULL;

   return;
}

//reads a graph from an input file
//The format of the input file is:
//number of vertices, followed by each edge
//It returns the graph as an incidence matrix
int** readgraph(char* filename,int& nvertices)
{
   FILE* in = fopen(filename,"r");
   if(NULL==in)
   {
      fprintf(stderr,"Cannot open input file [%s]\n",filename);
      exit(1);
   }

   int i,j;
   fscanf(in,"%d",&i);
   nvertices = i;

   int** graph = allocgraph(nvertices);

   while(2==fscanf(in,"%d %d",&i,&j))
   {
      graph[i-1][j-1] = graph[j-1][i-1] = 1;
   }

   fclose(in);

   return(graph);
}

//prints out a graph in Graphviz format
void printgraph(char* filename,int** graph,int nvertices)
{
   int i,j;

   FILE* out = fopen(filename,"w");
   if(NULL==out)
   {
      fprintf(stderr,"Cannot open output file [%s]\n",filename);
      exit(1);
   }

   fprintf(out,"graph G {\n");
   for(i=0;i<nvertices-1;i++)
   {
      for(j=i+1;j<nvertices;j++)
      {
	 if(graph[i][j])
	 {
	    fprintf(out,"%d -- %d;\n",i+1,j+1);
         }
      }
   }
 
   fprintf(out,"}\n");
   fclose(out);
   return;
}

void findneighbor(int** graph,int* record, int nvertices, int myvertex)
{
    int i;
	for(i=0;i<nvertices;i++){
	     if(graph[myvertex][i]==1)
	    {
		    if(record[i]==0){
			record[i]=1;
			findneighbor(graph,record,nvertices,i);
			}
			else continue;
		
		}
	}
return;
}

void findConComp(int myvertex,int** graph,int nvertices){
     int* record=new int[nvertices];
	 int i;
	 for(i=0;i<nvertices;i++){
	     record[i]=0;
	 }
	 record[myvertex]=1;
	 #initialize myvertex record as 1 to avoided rescursive 
	 
	 #find possible connected components of myvertex;
     findneighbor(graph,record,nvertices,myvertex);
	 
	 #print all the connected components
	 printf("The connected components of vertex %d: are: ",myvertex+1);
	 for(i=0;i<nvertices;i++){
	 if((i!=myvertex)&&(record[i]==1)) printf(" %d",i+1);
	 }
	 printf("\n");
	 
	 #free vector pointer record
	 free(record);
     return;
}