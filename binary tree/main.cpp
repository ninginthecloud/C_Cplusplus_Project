/*
 This program sorts the elements of a given vector.
*/

#include "vectors.h"
#include "tree.h"

int main()
{
	int i;
	int n = 8;
	
	char InputFile[] = "somenumbers.txt";
     char OutputFile[] = "mybinarytree.dot";

	//initialise a vector of length n
	double* vector = allocvector(n);
	
	readvector(n,vector,InputFile);
 
        printf("Original vector\n");
        printvector(vector,n);
        printf("\n\n");

        //create a binary tree
        LPNode mytree = MakeNewNode(vector[0]);
        //now add all the other elements of the vector
        for(i=1;i<n;i++)
	{
	  treeInsert(mytree,vector[i]);
        }

	//sort in increasing order by traversing the tree in inorder
	double* vectorIncreasing = allocvector(n);
    

	printf("Put %d number in increasing order:\n",n);
   for(i=0;i<n;i++)
	{
	    printf("-------------------------------------------\n");
        LPNode MinNode=FindMinValue(mytree);//find minnode
        vectorIncreasing[i]=MinNode->key;
        printf("   %4.2f\n",vectorIncreasing[i]);
        mytree = DeleteMinNode(mytree,MinNode);	//delete minnode
		}

	//free the memory
	freevector(vector);
	freevector(vectorIncreasing);

	return(1);
}

