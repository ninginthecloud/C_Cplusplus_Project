#include "matrices.h"

//allocates the memory for a matrix with 
//n rows and p columns
double ** allocmatrix(int n,int p)
{
	int i;
	double** m;
	
	m = new double*[n];
	for(i=0;i<n;i++)
	{
		m[i] = new double[p];
		memset(m[i],0,p*sizeof(double));
	}
	return(m);
}

//frees the memory for a matrix with n rows
void freematrix(int n,double** m)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		delete[] m[i]; m[i] = NULL;
	}
	delete[] m; m = NULL;
	return;
}

//creates the copy of a matrix with n rows and p columns
void copymatrix(int n,int p,double** source,double** dest)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			dest[i][j] = source[i][j];
		}
	}
	return;
}

//reads from a file a matrix with n rows and p columns
void readmatrix(char* filename,int n,int p,double* m[])
{
	int i,j;
	double s;
	FILE* in = fopen(filename,"r");
	
	if(NULL==in)
	{
		printf("Cannot open input file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			fscanf(in,"%lf",&s);
			m[i][j] = s;
		}
	}
	fclose(in);
	return;
}

//prints the elements of a matrix in a file
void printmatrix(char* filename,int n,int p,double** m)
{
	int i,j;
	double s;
	FILE* out = fopen(filename,"w");
	
	if(NULL==out)
	{
		printf("Cannot open output file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		fprintf(out,"%.3lf",m[i][0]);
		for(j=1;j<p;j++)
		{
			fprintf(out,"\t%.3lf",m[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	return;
}

//creates the transpose of the matrix m
double** transposematrix(int n,int p,double** m)
{
	int i,j;
	
	double** tm = allocmatrix(p,n);
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<n;j++)
		{
			tm[i][j] = m[j][i];
		}
	}	
	
	return(tm);
}

//calculates the dot (element by element) product of two matrices m1 and m2
//with n rows and p columns; the result is saved in m
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			m[i][j] = m1[i][j]*m2[i][j];
		}
	}
	
	return;
}

//calculates the product of a nxp matrix m1 with a pxl matrix m2
//returns a nxl matrix m
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m)
{
	int i,j,k;
	double s;
	
	for(i=0;i<n;i++)
	{
		for(k=0;k<l;k++)
		{
			s = 0;
			for(j=0;j<p;j++)
			{
				s += m1[i][j]*m2[j][k];
			}
			m[i][k] = s;
		}
	}
	return;
}

void set_mat_identity(int p, double *A)
{
 int i;

 for(i = 0; i < p * p; i++) A[i] = 0;
 for(i = 0; i < p; i++) A[i * p + i] = 1;
 return;
}

//computes the inverse of a symmetric positive definite matrix
void inverse(int p,double** m)
{
  int i,j,k;
  double* m_copy = (double*)malloc((p * p) * sizeof(double));
  double* m_inv = (double*)malloc((p * p) * sizeof(double));

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m_copy[k] = m[i][j];
        k++;
     }
  }

  set_mat_identity(p, m_inv);

  //-----  Use LAPACK  -------
  if(0!=(k=clapack_dposv(CblasRowMajor, CblasUpper, p, p, m_copy, p, m_inv, p)))
  {
    fprintf(stderr,"Something was wrong with clapack_dposv [%d]\n",k);
     exit(1);
  }
  //--------------------------

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m[i][j] = m_inv[k];
        k++;
     }
  }  

  free(m_copy);
  free(m_inv);

  return;
}


//computes the log of the determinant of a symmetric positive definite matrix
double logdet(int p,double** m)
{
	int i,j;
	char jobvl = 'N';
	char jobvr = 'N';
	int lda = p;
	double wr[2*p];
	double wi[2*p];
	double vl[p][p];
	int ldvl = p*p;
	double vr[p][p];
	int ldvr = p*p;
	double work[p*p];
	int lwork = p*p;
	double a[p][p];
	int info;
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			a[i][j] = m[i][j];
		}
	}
	dgeev_(&jobvl,&jobvr,&p,(double*)a,&lda,(double*)wr,(double*)wi,(double*)vl, 
		  &ldvl,(double*)vr,&ldvr,(double*)work,&lwork,&info);

	if(0!=info)
	{
		printf("Smth wrong in the call of 'dgeev' error is [info = %d]\n",info);
		exit(1);
	}	   
	
	double logdet = 0;
	for(i=0;i<p;i++) logdet+=log(wr[i]);	
	return(logdet);
}

//marglik
double marglik(int n,int p,double** data,int lenA,int* A)
{


//first 
//define constant value
 double const_value1;
 const_value1=lgamma((2.0+n+lenA)/2);
 double const_value2;
 const_value2=lgamma((lenA+2.0)/2);
   
 //second
 //define Da, which store the covariates values;
double** Da=allocmatrix(n,lenA);
int row;
int col;
int index;
  for(row=0;row<n;row++)
  {
    for(col=0;col<lenA;col++)
    {
	 index=A[col]-1;
     Da[row][col]=data[row][index];
     }
  } 
 
 //third 
 //define lenA dim identity matrix;

 double** I_lenA=allocmatrix(lenA,lenA);

	for (row=0;row<lenA;row++)
	{
		for(col=0;col<lenA;col++)
		{
			if(col==row){I_lenA[row][col]=1;}
			else{I_lenA[row][col]=0;}
		}
	}


//fourth
//define transpose matrix Da;
double** TDa=allocmatrix(lenA,n);
TDa=transposematrix(n,lenA,Da);

//fifth
//define ProdDa
double** ProdDa=allocmatrix(lenA,lenA);
matrixproduct(lenA,n,lenA,TDa,Da,ProdDa);
//sixth
//define Ma, which is the sum of I_lenA and ProdDa; 
double** Ma=allocmatrix(lenA,lenA);
 for(row=0;row<lenA;row++)
  {
    for(col=0;col<lenA;col++)
    {
     Ma[row][col]=I_lenA[row][col]+ProdDa[row][col];
     }
  } 

  //seventh
  //define log determinant of Ma
 double logdetMa;
 logdetMa=logdet(lenA,Ma);
 
 //eighth
 //define D_1, which is the response vector
 double** D1=allocmatrix(n,1);
 copymatrix(n,1,data,D1);
/*for (row=0; row<n; row++)
	{
		D1[row][0]=data[row][0];
	}
*/
//ninth
//define transpose matrix D_1
 double** TD1=allocmatrix(1,n);
 TD1=transposematrix(n,1,D1);
//tenth
//define inverse Ma
double** InverseMa=allocmatrix(lenA,lenA);
 copymatrix(lenA,lenA,Ma,InverseMa);//do I need to copy matrix, since I have allocate memory?
 inverse(lenA,InverseMa);

 //eleventh
 //define TDaD1;
 double** TDaD1=allocmatrix(lenA,1);
 matrixproduct(lenA,n,1,TDa,D1,TDaD1);
 
 //twelfthth
 //define TD1Da;
 double** TD1Da=allocmatrix(1,lenA);
 TD1Da=transposematrix(lenA,1,TDaD1);
 
 //thirteenth
 //define product TD1D1;
 double** TD1D1=allocmatrix(1,1);
 matrixproduct(1,n,1,TD1,D1,TD1D1);
 
 //fourteenth
 //define product matrix TD1Da and InverseMa;
 double** prod_right=allocmatrix(1,lenA);
 matrixproduct(1,lenA,lenA,TD1Da,InverseMa,prod_right);
 
 //fifteenth
 //define product value prod_right and TDaD1;
 double** prod_final=allocmatrix(1,1);
 matrixproduct(1,lenA,1,prod_right,TDaD1,prod_final);
 
 //sixteenth
 //define constant value const_value3
 double const_value3;
 const_value3=1+TD1D1[0][0]-prod_final[0][0];
 //
 /* for debug
    printmatrix("Da.txt",n,lenA,Da);
	printmatrix("I_lenA.txt",lenA,lenA,I_lenA);
	printmatrix("TDa.txt",lenA,n,TDa);
	printmatrix("ProdDa.txt",lenA,lenA,ProdDa);
	printmatrix("Ma.txt",lenA,lenA,Ma);
	printmatrix("D1.txt",n,1,D1);
    printmatrix("TD1.txt",1,n,TD1);
	printmatrix("InverseMa.txt",lenA,lenA,InverseMa);
	printmatrix("TDaD1.txt",lenA,1,TDaD1);
	printmatrix("TD1Da.txt",1,lenA,TD1Da);
	printmatrix("TD1D1.txt",1,1,TD1D1);
	printmatrix("prod_right.txt",1,lenA,prod_right);
	printmatrix("prod_final.txt",1,1,prod_final);
	*/	
 //freematrix
 //free matrix
freematrix(n,Da);//second
freematrix(lenA,I_lenA);//third
freematrix(lenA,TDa);//fourth
freematrix(lenA,ProdDa);//fifth
freematrix(lenA,Ma);//sixth
freematrix(n,D1);//eighth
freematrix(1,TD1);//ninth
freematrix(lenA,InverseMa);//tenth
freematrix(lenA,TDaD1);//eleventh
freematrix(1,TD1Da);//twelfth
freematrix(1,TD1D1);//thirteenth
freematrix(1,prod_right);//fourteenth
freematrix(1,prod_final);//fifteenth


 //seventeenth
 //define final return value;
 double result;
 result= const_value1-const_value2-0.5*logdetMa-(n+lenA+2.0)/2*log(const_value3);
 return(result);


}

