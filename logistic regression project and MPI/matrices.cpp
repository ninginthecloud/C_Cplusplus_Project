#include "matrices.h"



//creates the transpose of the matrix m
gsl_matrix* transposematrix(gsl_matrix* m)
{
	int i,j;
	
	gsl_matrix* tm = gsl_matrix_alloc(m->size2,m->size1);
	
	for(i=0;i<tm->size1;i++)
	{
		for(j=0;j<tm->size2;j++)
		{
		  gsl_matrix_set(tm,i,j,gsl_matrix_get(m,j,i));
		}
	}	
	
	return(tm);
}

//calculates the product of a nxp matrix m1 with a pxl matrix m2
//returns a nxl matrix m
void matrixproduct(gsl_matrix* m1,gsl_matrix* m2,gsl_matrix* m)
{
	int i,j,k;
	double s;
	
	for(i=0;i<m->size1;i++)
	{
	  for(k=0;k<m->size2;k++)
	  {
	    s = 0;
	    for(j=0;j<m1->size2;j++)
	    {
	      s += gsl_matrix_get(m1,i,j)*gsl_matrix_get(m2,j,k);
	    }
	    gsl_matrix_set(m,i,k,s);
	  }
	}
	return;
}


//computes the inverse of a positive definite matrix
//the function returns a new matrix which contains the inverse
//the matrix that gets inverted is not modified
gsl_matrix* inverse(gsl_matrix* K)
{
	int j;
	
	gsl_matrix* copyK = gsl_matrix_alloc(K->size1,K->size1);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(copyK,K))
	{
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}
	
	gsl_matrix* inverse = gsl_matrix_alloc(K->size1,K->size1);
	gsl_permutation *myperm = gsl_permutation_alloc(K->size1);
	
	if(GSL_SUCCESS!=gsl_linalg_LU_decomp(copyK,myperm,&j))
	{
		printf("GSL failed LU decomposition.\n");
		exit(1);
	}
	if(GSL_SUCCESS!=gsl_linalg_LU_invert(copyK,myperm,inverse))
	{
		printf("GSL failed matrix inversion.\n");
		exit(1);
	}
	gsl_permutation_free(myperm);
	gsl_matrix_free(copyK);
	
	return(inverse);
}


//FUNCTION DEFINITION
//Newton-Raphson algothrim
//Loglikelihood 

//Logdet 
//computes the log of the determinant of a symmetric positive definite matrix
double logdet(gsl_matrix* K)
{
        int i;

	gsl_matrix* CopyOfK = gsl_matrix_alloc(K->size1,K->size2);
	gsl_matrix_memcpy(CopyOfK,K);
	gsl_permutation *myperm = gsl_permutation_alloc(K->size1);
	if(GSL_SUCCESS!=gsl_linalg_LU_decomp(CopyOfK,myperm,&i))
	{
		printf("GSL failed LU decomposition.\n");
		exit(1);
	}
	double logdet = gsl_linalg_LU_lndet(CopyOfK);
	gsl_permutation_free(myperm);
	gsl_matrix_free(CopyOfK);
	return(logdet);
}
//inverse of the logit function
double inverseLogit(double x)
{
//double result=gsl_sf_exp((-1.0)*gsl_sf_log_1plusx(gsl_sf_exp((-1)*x)));
double result=1.0/(1.0+exp((-1.0)*x));
return(result);
}

//function for the computation of the hessian
double inverseLogit2(double x)
{
double result= 1.0/(1.0+exp((-1.0)*x))*1.0/(1.0+exp(x));
return(result);
}
 
 //computes pi_i = P(y_i = 1 | x_i)
 gsl_matrix* getPi(gsl_vector* x, gsl_matrix* beta)
 {
  //X0=(1,x), product X0*beta beta is a 2-dim vector
  int i;
  int p=x->size;
  gsl_matrix* A=gsl_matrix_calloc(p,2);
  gsl_matrix* y=gsl_matrix_calloc(p,1);
  gsl_matrix* Pi=gsl_matrix_calloc(p,1);
  for(i=0;i<p;i++){
     gsl_matrix_set(A,i,0,1.0);
	 gsl_matrix_set(A,i,1,gsl_vector_get(x,i));
  }
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0,A,beta,0.0,y);
    for(i=0;i<p;i++)
  {
     gsl_matrix_set(Pi,i,0,inverseLogit(gsl_matrix_get(y,i,0)));
  }
  gsl_matrix_free(A);
  gsl_matrix_free(y);
  return(Pi);
   }
 
//another function for the computation of the Hessian
 gsl_matrix* getPi2(gsl_vector* x, gsl_matrix* beta)
 {
  //X0=(1,x), product X0*beta beta is a 2-dim vector
  int i;
  int p=x->size;
  gsl_matrix* A=gsl_matrix_calloc(p,2);
  gsl_matrix* y=gsl_matrix_calloc(p,1);
  gsl_matrix* Pi=gsl_matrix_calloc(p,1);
  for(i=0;i<p;i++){
     gsl_matrix_set(A,i,0,1.0);
	 gsl_matrix_set(A,i,1,gsl_vector_get(x,i));
  }
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0,A,beta,0.0,y);
    for(i=0;i<p;i++)
  {
     gsl_matrix_set(Pi,i,0,inverseLogit2(gsl_matrix_get(y,i,0)));
  }
  gsl_matrix_free(A);
  gsl_matrix_free(y);
  return(Pi);
   }
 
//logistic log-likelihood (formula (3) in your handout)
//gsl_sf_log_1plusx compute \log(1 + x) for x > -1 using an algorithm that is accurate for small x. 
double logisticLoglik(gsl_vector* y,gsl_vector* x,gsl_matrix* beta)
{
  gsl_matrix* Pi = getPi(x,beta);
  double logisticL=0.0;
  int i;
  for(i=0;i<y->size;i++){
    double pi=gsl_matrix_get(Pi,i,0);
	double yi=gsl_vector_get(y,i);
    logisticL=logisticL+yi*gsl_sf_log(pi)+(1.0-yi)*gsl_sf_log(1.0-pi);
  }
  gsl_matrix_free(Pi);
  return(logisticL);
} 

//calculates l^*(\beta_0,\beta_1)//!pi?
double lStar(gsl_vector* y,gsl_vector* x,gsl_matrix* beta)
{
  double prodbeta=pow((gsl_matrix_get(beta,0,0)),2.0)+pow((gsl_matrix_get(beta,1,0)),2.0);
  double result=logisticLoglik(y,x,beta)-0.5*prodbeta-gsl_sf_log(M_PI*2);
  return(result);
}
//obtain the gradient for Newton-Raphson
//gsl_blas_daxpy is used to calculate vector y+Pi
//
gsl_matrix* getGradient(gsl_vector* y,gsl_vector* x,gsl_matrix* beta)
{
  gsl_matrix* gradient= gsl_matrix_calloc(beta->size1,1);
  gsl_matrix* Pi = getPi(x,beta); 
  //generate design matrix x0
  gsl_matrix* x0 = gsl_matrix_calloc(x->size,2);
  gsl_matrix_set_all(x0,1.0);
  gsl_matrix_set_col(x0,1,x);
  // allocate a temp matrix to store value
  gsl_matrix* temp = gsl_matrix_calloc(Pi->size1,1);
  gsl_matrix_set_col(temp,0,y);
  gsl_matrix_sub(temp,Pi);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,x0,temp,0.0,gradient);
  gsl_matrix_sub(gradient,beta);
  //free memory
  gsl_matrix_free(Pi);
  gsl_matrix_free(x0);
  gsl_matrix_free(temp);
  
  return(gradient);
}
 
//obtain the Hessian for Newton-Raphson
gsl_matrix* getHessian(gsl_vector* y, gsl_vector* x,gsl_matrix* beta)
{
  gsl_matrix* hessian = gsl_matrix_calloc(beta->size1,beta->size1);
  gsl_matrix* Pi2 = getPi2(x,beta);
  //generate design matrix x0
  gsl_matrix* x0 = gsl_matrix_calloc(x->size,2);
  gsl_matrix_set_all(x0,1.0);
  gsl_matrix_set_col(x0,1,x);
  //generate identity matrix
  gsl_matrix* identity = gsl_matrix_calloc(2,2);
  gsl_matrix_set_all(identity,1.0);
  gsl_matrix_set_identity(identity);
  gsl_vector_const_view Pi2view = gsl_matrix_const_column(Pi2,0);
  gsl_matrix* temp = gsl_matrix_alloc(x->size,2);
  gsl_matrix_set_col(temp,0,&Pi2view.vector);
  gsl_matrix_set_col(temp,1,&Pi2view.vector);
  gsl_matrix_mul_elements(temp,x0);
  //crossprod
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,temp,x0,0.0,hessian);
  gsl_matrix_sub(hessian,identity);
  gsl_matrix_scale(hessian,-1.0);
  //memory free
  gsl_matrix_free(Pi2);
  gsl_matrix_free(x0);
  gsl_matrix_free(identity);
  gsl_matrix_free(temp);
 
  return(hessian);
}
  //this function implements our own Newton-Raphson procedure
gsl_matrix* getcoefNR(int response,int explanatory,gsl_matrix* data)
{
  double epsilon=1e-6;
  //2-dim vector of coefficients`
  gsl_matrix* beta = gsl_matrix_calloc(2,1);//here need to be improved, dont manually allocate
  gsl_vector_view y = gsl_matrix_column(data,(response-1));
  gsl_vector_view x = gsl_matrix_column(data,(explanatory-1));
  //current value of log-likelihood
  double currentLStar = lStar(&y.vector,&x.vector,beta);
  //infinite loop unless we stop it someplace inside
  int Iteration=0;
  while(1)
  {
    Iteration=Iteration+1;
    gsl_matrix* Hessian=getHessian(&y.vector,&x.vector,beta);
	gsl_matrix* inverseHessian=inverse(Hessian);
	gsl_matrix* Gradiant=getGradient(&y.vector,&x.vector,beta);
	gsl_matrix* a=gsl_matrix_calloc(2,1);
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, inverseHessian,Gradiant,0.0,a);
	gsl_matrix* newbeta=gsl_matrix_calloc(beta->size1,beta->size2);
	gsl_matrix_memcpy(newbeta,beta);
    gsl_matrix_sub(newbeta,a);  
    // get new loglikelihood 
    double newLStar	= lStar(&y.vector,&x.vector,newbeta);
    
    //at each iteration the log-likelihood must increase
    if(newLStar<currentLStar)
    {
      printf("CODING ERROR!! %d\t  %.51f\t %.51f\n",Iteration,newLStar,currentLStar);
      //free memory;
      gsl_matrix_free(Hessian);
      gsl_matrix_free(inverseHessian);
      gsl_matrix_free(Gradiant);
      gsl_matrix_free(a);
      gsl_matrix_free(newbeta);
     
	  break;
    }
    //stop if the log-likelihood does not improve by too much
    
	gsl_matrix* diff = gsl_matrix_calloc(beta->size1,beta->size2);
	gsl_matrix_memcpy(diff,newbeta);
	gsl_matrix_sub(diff,beta);
	
	int i;
	for(i=0;i<diff->size1;i++)
	{
	    gsl_matrix_set(diff,i,0,fabs(gsl_matrix_get(diff,i,0)));
	}
	
	if(gsl_matrix_max(diff)<epsilon)
    {
	//free memory
	//free memory;
  gsl_matrix_free(Hessian);
  gsl_matrix_free(inverseHessian);
  gsl_matrix_free(Gradiant);
  gsl_matrix_free(a);
  gsl_matrix_free(newbeta);
  gsl_matrix_free(diff);
      break; 
    }
	gsl_matrix_memcpy(beta,newbeta);
    currentLStar = newLStar;
	//free memory;
  gsl_matrix_free(Hessian);
  gsl_matrix_free(inverseHessian);
  gsl_matrix_free(Gradiant);
  gsl_matrix_free(a);
  gsl_matrix_free(newbeta);
  gsl_matrix_free(diff);  }
  
  return(beta);
}

//creates the Cholesky decomposition of a matrix
gsl_matrix* makeCholesky(gsl_matrix* K)
{
	int i,j;
	
	gsl_matrix* Phi = gsl_matrix_alloc(K->size1,K->size1);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(Phi,K))
	{
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}
	if(GSL_SUCCESS!=gsl_linalg_cholesky_decomp(Phi))
	{
		printf("GSL failed Cholesky decomposition.\n");
		exit(1);
	}
	for(i=0;i<Phi->size1;i++)
	{
		for(j=i+1;j<Phi->size2;j++)
		{
			gsl_matrix_set(Phi,i,j,0.0);
		}
	}
	
	return(Phi);
}

//samples from the multivariate normal distribution N(Mu,Sigma)
//the samples are saved in the matrix "Samples"
void randomMVN(gsl_rng* r,gsl_matrix* Samples,gsl_matrix* Mu,gsl_matrix* Sigma)
{
  int i;
  gsl_matrix* Psi = makeCholesky(Sigma);
  gsl_matrix* Z = gsl_matrix_alloc(Sigma->size1,1);
  gsl_matrix* X = gsl_matrix_alloc(Sigma->size1,1);
  
  for(int asample=0;asample<Samples->size2;asample++)
  {
    for(i=0;i<Sigma->size1;i++)
    {
      gsl_matrix_set(Z,i,0,gsl_ran_ugaussian(r));
    }
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		    1.0, Psi, Z,
		    0.0, X); 
			
    //record the sample we just generated
    for(i=0;i<Sigma->size1;i++)
    {
      gsl_matrix_set(Samples,i,asample,
		     (gsl_matrix_get(X,i,0)+gsl_matrix_get(Mu,i,0)));
    }
  }
  
  //free memory
  gsl_matrix_free(Psi);
  gsl_matrix_free(Z);
  gsl_matrix_free(X);
  return;
}

//performs one iteration of the Metropolis-Hastings algorithm
void mhLogisticRegression(gsl_rng* r,gsl_vector* y,gsl_vector* x,
                                   gsl_matrix* beta,gsl_matrix* invNegHessian)
{
  
  gsl_matrix* betaCandidate = gsl_matrix_calloc(beta->size1,beta->size2);
  randomMVN(r,betaCandidate,beta,invNegHessian);
  
  //current value of lstar
  double currentLStar = lStar(y,x,beta);

  //values of lstar associated with betaCandidate
  double candidateLStar = lStar(y,x,betaCandidate);

/* 
 if(candidateLStar>=currentLStar)
  { 
    return(betaCandidate); 
  }
  */
  double u = gsl_rng_uniform_pos(r);
  if(log(u)<=(candidateLStar-currentLStar))//
  {
   gsl_matrix_memcpy(beta,betaCandidate);
  }
  
  //reject the move and stay at the current state
  gsl_matrix_free(betaCandidate);
  return;
}

//Problem 1: approximate the marginal likelihood (3) using the Laplace approximation
//'betaMode' is the posterior mode obtained using Newton-Raphson
double getLaplaceApprox(int response,int explanatory,
                                     gsl_matrix* data,gsl_matrix* betaMode)
{
  gsl_vector_view y = gsl_matrix_column(data,(response-1));
  gsl_vector_view x = gsl_matrix_column(data,(explanatory-1));
  
  //current value of log-likelihood
  double maxLogLik = logisticLoglik(&y.vector,&x.vector,betaMode);
  double prodbeta=pow((gsl_matrix_get(betaMode,0,0)),2.0)+pow((gsl_matrix_get(betaMode,1,0)),2.0);
  gsl_matrix* resultgetHessian=getHessian(&y.vector,&x.vector,betaMode);
  gsl_matrix_scale(resultgetHessian,-1.0);
  double logmarglik = maxLogLik-0.5*logdet(resultgetHessian)-0.5*prodbeta;
  //free memory
  gsl_matrix_free(resultgetHessian);
  return(logmarglik);
}

//Numerically calculate the marginal likelihood (3) using Monte Carlo (numerically stable version)
double getMonteCarlo(gsl_rng* r,int response,int explanatory,gsl_matrix* data,int NumberOfIterationsMC)
{
   gsl_vector_view y = gsl_matrix_column(data,(response-1));
   gsl_vector_view x = gsl_matrix_column(data,(explanatory-1));

   gsl_vector* loglikVec = gsl_vector_calloc(NumberOfIterationsMC);
   int i;
  for(i=0;i<NumberOfIterationsMC;i++)
  {
    gsl_matrix* tempbeta=gsl_matrix_calloc(2,1);
    gsl_matrix_set(tempbeta,0,0,gsl_ran_ugaussian(r));
	gsl_matrix_set(tempbeta,1,0,gsl_ran_ugaussian(r));
    gsl_vector_set(loglikVec,i,logisticLoglik(&y.vector,&x.vector,tempbeta));
	gsl_matrix_free(tempbeta);
 }
  double maxloglikVec = gsl_vector_max(loglikVec);
  for(i=0;i<NumberOfIterationsMC;i++){
    double tempexp = gsl_vector_get(loglikVec,i)-maxloglikVec;
    gsl_vector_set(loglikVec,i,gsl_sf_exp(tempexp));
  }
  double tempmean=gsl_stats_mean(loglikVec->data,loglikVec->stride,loglikVec->size);
  double logmarglik = gsl_sf_log(tempmean)+maxloglikVec;
  //free memory
  gsl_vector_free(loglikVec);
  return(logmarglik);
}

//Problem 2 calculate the posterior means of beta0 and beta1 by
//sampling from the joint posterior distribution of beta0 and beta1
//with the Metropolis-Hastings algorithm
gsl_matrix* getPosteriorMeans(gsl_rng*r,int response,int explanatory,gsl_matrix* data,
                                   gsl_matrix* betaMode,int NumberOfIterations)
{
    gsl_vector_view y = gsl_matrix_column(data,(response-1));
    gsl_vector_view x = gsl_matrix_column(data,(explanatory-1));
    
	gsl_matrix* betaBayes = gsl_matrix_calloc(betaMode->size1,1);
    gsl_matrix* betaCurrent = gsl_matrix_calloc(betaMode->size1,1);
    gsl_matrix_memcpy(betaCurrent,betaMode);
   
   //the negative of the inverse Hessian matrix evaluated at the posterior mode
   gsl_matrix* resultgetHessian=getHessian(&y.vector,&x.vector,betaMode);  
   gsl_matrix* invNegHessian=inverse(resultgetHessian);
   gsl_matrix_scale(invNegHessian,(-1.0));
   int iteration; 
   for(iteration=0;iteration<NumberOfIterations;iteration++)
   {
      mhLogisticRegression(r,&y.vector,&x.vector,betaCurrent,invNegHessian);
   	  gsl_matrix_add(betaBayes,betaCurrent);
	  }
   gsl_matrix_scale(betaBayes,1.0/NumberOfIterations);
   //free memory
   gsl_matrix_free(betaCurrent);
   gsl_matrix_free(resultgetHessian);
   gsl_matrix_free(invNegHessian);
   return(betaBayes);
}
 
//bayesLogistic function
gsl_vector* bayesLogistic(gsl_rng* r,int predictor,int response,gsl_matrix* data,
             int NumberOfIterations,int NumberOfIterationsMC)
{
    //use Newton-Raphson to calculate the mode of the posterior distribution (2)
    //the mode is needed in the Laplace approximation and in the Metropolis-Hastings algorithm
    gsl_matrix* tempbetaMode = getcoefNR(response,predictor,data);

    //Problem 1: approximate the marginal likelihood (3) using the Laplace approximation
   // double templogmarglik = getLaplaceApprox(response,predictor,data,tempbetaMode);
  
    //Numerically calculate the marginal likelihood (3) using Monte Carlo (stable version)
    double templogmarglikMC  = getMonteCarlo(r,response,predictor,data,NumberOfIterationsMC);

    //Problem 2: calculate the posterior means of beta0 and beta1
    //by sampling from the joint posterior of beta0 and beta1
    //using the Metropolis-Hastings algorithm
    gsl_matrix* tempbetaBayes = getPosteriorMeans(r,response,predictor,data,tempbetaMode,NumberOfIterations);
   
    //Problem 1: approximate the marginal likelihood (3) using the Laplace approximation
    double templogmarglik = getLaplaceApprox(response,predictor,data,tempbetaBayes);
	//create new element 
	gsl_vector* result=gsl_vector_calloc(5);
	
	gsl_vector_set(result,0,(double)predictor);
	gsl_vector_set(result,1,templogmarglik);
	gsl_vector_set(result,2,templogmarglikMC);
	gsl_vector_set(result,3,gsl_matrix_get(tempbetaBayes,0,0));
	gsl_vector_set(result,4,gsl_matrix_get(tempbetaBayes,1,0));
		
	//free memory
	gsl_matrix_free(tempbetaBayes);
	gsl_matrix_free(tempbetaMode);
	return(result);
   }

 
 
 
 
 