#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

//THESE ARE GSL FUNCTIONS
//YOU DO NOT NEED TO INCLUDE ALL THESE HEADER FILES IN YOUR CODE
//JUST THE ONES YOU ACTUALLY NEED;
//IN THIS APPLICATION, WE ONLY NEED gsl/gsl_matrix.h
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>


gsl_matrix* transposematrix(gsl_matrix* m);
void matrixproduct(gsl_matrix* m1,gsl_matrix* m2,gsl_matrix* m);
gsl_matrix* inverse(gsl_matrix* K);
double logdet(gsl_matrix* K);
double inverseLogit(double x);
double inverseLogit2(double x);
gsl_matrix* getPi(gsl_vector* x, gsl_matrix* beta);
gsl_matrix* getPi2(gsl_vector* x, gsl_matrix* beta);
double logisticLoglik(gsl_vector* y,gsl_vector* x,gsl_matrix* beta);
double lStar(gsl_vector* y,gsl_vector* x,gsl_matrix* beta);
gsl_matrix* getGradient(gsl_vector* y,gsl_vector* x,gsl_matrix* beta);
gsl_matrix* getHessian(gsl_vector* y, gsl_vector* x,gsl_matrix* beta);
gsl_matrix* getcoefNR(int response,int explanatory,gsl_matrix* data);
gsl_matrix* makeCholesky(gsl_matrix* K);
void randomMVN(gsl_rng* r,gsl_matrix* Samples,gsl_matrix* Mu,gsl_matrix* Sigma);
void mhLogisticRegression(gsl_rng* r,gsl_vector* y,gsl_vector* x,gsl_matrix* beta,gsl_matrix* invNegHessian);
double getLaplaceApprox(int response,int explanatory,
                                     gsl_matrix* data,gsl_matrix* betaMode);
double getMonteCarlo(gsl_rng* r,int response,int explanatory,gsl_matrix* data,int NumberOfIterations);
gsl_matrix* getPosteriorMeans(gsl_rng* r,int response,int explanatory,gsl_matrix* data,
                                   gsl_matrix* betaMode,int NumberOfIterations);
gsl_vector* bayesLogistic(gsl_rng* r,int predictor,int response,gsl_matrix* data,
             int NumberOfIterations,int NumberOfIterationsMC);
