/*
 This program computes the discrete probability of 0-(n) Compile the program using the makefile provided.
 
 Run the program using the command:

 mpirun -np 10 parallel
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include "matrices.h"
#include "Bayes.h"

/////////////////////////////////////////////////////
//MPI MESSAGES TO THE SLAVES                       //
//the master process can ask the slaves to estimate//
//another pi=P(X=i) or to die.                     //
/////////////////////////////////////////////////////
#define GETPI 1
#define DIETAG 0


void master(int n)
{
  int nMaxReg=5;
  LPBayes reg=new Bayes;
  reg->Next=NULL;
  
  int i;
  int rank;
  int ntasks;
  int jobsRunning;
  char outfilename[] = "best5.txt";
  int work[1]; //used to transmit a work request to the slaves
  double workresults[5]; //used to receive results from the slaves
  MPI_Status status;	// MPI information

  // Find out how many slaves there are
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  fprintf(stderr, "Total Number of processors = %d\n",ntasks);

  jobsRunning = 1;
  for(i=1;i<=n;i++)
  {
    work[0] = i; //work request to estimate pi=P(X=i)
    if(jobsRunning < ntasks) // Do we have an available processor?
    {
      // Send out a work request
      MPI_Send(&work, 	// the vector with the variable
	       1, 	// the size of the vector
	       MPI_INT,	// the type of the vector
               jobsRunning,	// the ID of the slave to use
               GETPI,	// tells the slave what to do
               MPI_COMM_WORLD); 

       printf("Master sends out work request [%d] to slave [%d]\n",
              work[0],jobsRunning);

       // Increase the # of processors in use
       jobsRunning++;
      }
      else // all the processors are in use!
      {
	 //RECEIVE RESULTS
         MPI_Recv(workresults,	// where to store the results
 		  5,		// the size of the vector
		  MPI_DOUBLE,	// the type of the vector
	 	  MPI_ANY_SOURCE,
		  MPI_ANY_TAG, 	
		  MPI_COMM_WORLD,
		  &status);
	 printf("Master has received the result of work request [%d] from slave [%d]\n",
                (int)workresults[0],status.MPI_SOURCE);

	 //STORE RESULTS in a link!!!TODO
	 AddRegression(nMaxReg,reg,(int)workresults[0],
                   workresults[1],workresults[2],workresults[3],workresults[4]);
     

	 printf("Master sends out work request [%d] to slave [%d]\n",
                work[0],status.MPI_SOURCE);
	 // Send out a new work order to the processors that just
         // returned
         MPI_Send(&work,
                  1,
                  MPI_INT,
                  status.MPI_SOURCE, // the slave that just returned
                  GETPI,
                  MPI_COMM_WORLD);
      }
  }

  /////////////////////////////////////////////////////////
  //COLLECT RESULTS FOR ALL THE OUTSTANDING WORK REQUESTS//
  /////////////////////////////////////////////////////////
  for(rank=1; rank<jobsRunning; rank++)
  {
    //RECEIVE RESULTS
    MPI_Recv(workresults,	// where to store the results
 	     5,		// the size of the vector
	     MPI_DOUBLE,	// the type of the vector
	     MPI_ANY_SOURCE,
	     MPI_ANY_TAG, 	
	     MPI_COMM_WORLD,
	     &status);
     printf("Master has received the result of work request [%d] from slave [%d]\n",
            (int)workresults[0],status.MPI_SOURCE);

     //STORE RESULTS in a link list! TODO
	  AddRegression(nMaxReg,reg,(int)workresults[0],
                   workresults[1],workresults[2],workresults[3],workresults[4]);

     
  }

  printf("Master tells the slave to die\n");
  // Shut down the slave processes
  for(rank=1; rank<ntasks; rank++)
  {
    printf("Master is killing slave [%d]\n",rank);
    MPI_Send(0,
	     0,
             MPI_INT,
             rank,		// shutdown this particular node
             DIETAG,		// tell it to shutdown
	     MPI_COMM_WORLD);
  }

  ////////////////////
  //SHOW THE RESULTS//
  ////////////////////
   SaveRegressions(outfilename,reg);
   
  //free memory
  DeleteAllRegressions(reg);
  delete reg;
  
  return;
}

void slave(int NumberOfIterationsMC,int NumberOfIterations,int slavename)
{
  // read datafile
    int size = 148; //sample size
	int p = 61; //number of variables
	int response = 61;//response value comes from 61 column in data
	char datafilename[] = "534finalprojectdata.txt";
	//allocate the data matrix
	gsl_matrix* data = gsl_matrix_alloc(size,p);
    //read the data
	FILE* datafile = fopen(datafilename,"r");
	if(NULL==datafile)
	{
		fprintf(stderr,"Cannot open data file [%s]\n",datafilename);
		return;
	}
	if(0!=gsl_matrix_fscanf(datafile,data))
	{
		fprintf(stderr,"File [%s] does not have the required format.\n",datafilename);
		return;
	}
	fclose(datafile);
  
  int i;
  int work[1]; //used to receive a work request from the Master
  double workresults[5]; //used to transmit results back to the Master
  MPI_Status status; //used for MPI communication

	
  //initialize the GSL RNG
  const gsl_rng_type* T;
  gsl_rng* r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  //VERY IMPORTANT: MAKE SURE THE SLAVES DO NOT GENERATE
  //THE SAME SEQUENCE OF RANDOM NUMBERS
  //each slave generates a different number of U(0,1)
  for(i=0;i<slavename*NumberOfIterations;i++)
  {
    gsl_rng_uniform(r);
  }

  //////////////////////////////////////////////////////
  //THE SLAVE LISTENS FOR INSTRUCTIONS FROM THE MASTER//
  //////////////////////////////////////////////////////
  int notDone = 1;
  while(notDone)
  {
     printf("Slave %d is waiting\n",slavename);
     MPI_Recv(&work,		// the inputs from the master
	      1,		// the size of the inputs
	      MPI_INT,		// the type of the inputs
              0,		// from the MASTER node (rank=0)
              MPI_ANY_TAG,	// any type of order is fine
              MPI_COMM_WORLD,
              &status);
      printf("Slave %d just received smth\n",slavename);

      // switch on the type of work request
      switch(status.MPI_TAG)
      {
         case GETPI:
        {   printf("Slave %d has received work request [%d]\n",
                  slavename,work[0]);
          
	   //work!
	   gsl_vector* result = bayesLogistic(r,work[0],response,data,
             NumberOfIterations,NumberOfIterationsMC);
       
       
       workresults[0] = (double)work[0];
	   workresults[1] = gsl_vector_get(result,1);//log marginal likelihood
	   workresults[2] = gsl_vector_get(result,2);//MC log marginal likelihood
	   workresults[3] = gsl_vector_get(result,3);// Bayes estimator beta0
	   workresults[4] = gsl_vector_get(result,4);// Bayes estimator beta1
         printf("slave beta0 beta1 %f,%f\n",workresults[3],workresults[4]);
		//free memory
		gsl_vector_free(result);	   
         
           // Send the results
           MPI_Send(&workresults,
                    5,
                    MPI_DOUBLE,
                    0,		// send it to the master
                    0,		// doesn't need a TAG
                    MPI_COMM_WORLD);

            printf("Slave %d finished processing work request [%d]\n",
                   slavename,work[0]);
		
            break;
			
			}

         case DIETAG:
            printf("Slave %d was told to die\n",slavename);

         default:
            notDone = 0;
      }
  }

  //free memory
  gsl_rng_free(r);
  gsl_matrix_free(data);
  return;
}

////////////////////////////////////////////////////////////////////
//The program should be called by specifying the number of batches//
//and the number of iterations per batch in the command line.     //
//As such, argc should be 4.                                      //
//For example, to complete Problem 2, you make the following call://
// mpirun -np 6 randbin 10 25 1000                                //
//in order to simulate from the discrete distribution associated  //
//with the binomial coefficients for n = 10 in 25 batches of      //
//1000 samples each. You requested MPI to create 6 copies of your //
//program, hence you get 1 Master and 5 Slaves                    //
////////////////////////////////////////////////////////////////////
int main(int argc,char** argv)
{
        
	    int myrank;

	if(4!=argc)
	{
	  fprintf(stderr,"USAGE: mpirun -np <NumberOfProcesses> parallelBayes <n Number of predictors> <Number of Iterations for MC> <Number of Iterations for posterior calculation>\n");
	  return(0);
	}

	///////////////////////////
	// START THE MPI SESSION //
	///////////////////////////
	MPI_Init(&argc, &argv);

	/////////////////////////////////////
	// What is the ID for the process? //   
	/////////////////////////////////////
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

	int n = atoi(argv[1]);
	int NumberOfIterationsMC = atoi(argv[2]);
	int NumberOfIterations = atoi(argv[3]);

	if(myrank==0)
	{
	  fprintf(stderr,
		"You asked for n=%d, %d Iterations for MC and %d Iterations for posterior calculation.\n",
		n,NumberOfIterationsMC,NumberOfIterations);
	  master(n);
	}
	else
	{
	  slave(NumberOfIterationsMC,NumberOfIterations,myrank);
	}

	// Finalize the MPI session
	MPI_Finalize();

	return(1);
}

