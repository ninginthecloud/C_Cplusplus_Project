#ifndef _Bayes
#define _Bayes

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <gsl/gsl_vector.h>


typedef struct myBayes* LPBayes;
typedef struct myBayes Bayes;

struct myBayes
{
  int apredictor; // predictor added into bayeslogistic regression
  double logmarglik; // log marginal likelihood value_comp
  double logmarglikMC;
  gsl_vector* betaBayes;
  LPBayes Next; //link to the next bayeslogistic regression
};

void AddRegression(int nMaxReg,LPBayes regressions,int predictor,
                   double alogmarlik,double alogmarglikMC,double beta0,double beta1);
void DeleteAllRegressions(LPBayes regressions);
void DeleteLastRegression(LPBayes regressions);				   
void SaveRegressions(char* filename,LPBayes regressions);


#endif
