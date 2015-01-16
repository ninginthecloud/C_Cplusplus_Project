/*
 FILE: Bayes.CPP
*/

#include "Bayes.h"

void AddRegression(int nMaxRegs,LPBayes regressions,int predictor,
                   double alogmarglik,double alogmarglikMC,double beta0,double beta1)
{
  int i, j = 0;

  LPBayes p = regressions;
  LPBayes pnext = p->Next;

  while(NULL != pnext && j < nMaxRegs)
  {
     //go to the next element in the list if the current
     //regression has a larger log marginal likelihood than
     //the new regression A
     if(pnext->logmarglikMC > alogmarglikMC)
     {
        p = pnext;
        pnext = p->Next;
     }
     else //otherwise stop; this is where we insert the new regression
     {
        break;
     }
     j++;
  }

  // if we reached "nMaxRegs" we did not beat any of the top 10 with the new regression.
  // Otherwise we add it like normal.

  if(nMaxRegs == j)
  {
	  return;
  }

  //create a new element of the list
  LPBayes newp = new Bayes;
  newp->apredictor= predictor;
  newp->logmarglik = alogmarglik;
  newp->logmarglikMC = alogmarglikMC;
  newp->betaBayes = gsl_vector_calloc(2);
  gsl_vector_set(newp->betaBayes,0,beta0);
  gsl_vector_set(newp->betaBayes,1,beta1);
  //insert the new element in the list
  p->Next = newp;
  newp->Next = pnext;
   // now we move through the list until we either reach the end of it, or reach the
  // element just after the "nMaxRegs" element.
  while(j < nMaxRegs && NULL!=pnext)
  {
	  p = pnext;
	  pnext = p->Next;
	  j++;
  }
  // if we reach nMaxRegs, we have to discard the new worst element in the list.
  if(nMaxRegs == j)
  {
	  DeleteLastRegression(regressions);
  }

  return;
}

//this function deletes all the elements of the list
//with the head "regressions"
//remark that the head is not touched
void DeleteAllRegressions(LPBayes regressions)
{
  //this is the first regression
  LPBayes p = regressions->Next;
  LPBayes pnext;

  while(NULL!=p)
  {
    //save the link to the next element of p
    pnext = p->Next;

    //delete the element specified by p
    //first free the memory of the vector of regressors
    gsl_vector_free(p->betaBayes);
    p->Next = NULL;
    delete p;

    //move to the next element
    p = pnext;
  }

  return;
}

//this function deletes the last element of the list
//with the head "regressions"
//again, the head is not touched
void DeleteLastRegression(LPBayes regressions)
{
  //this is the element before the first regression
  LPBayes pprev = regressions;
  //this is the first regression
  LPBayes p = regressions->Next;

  //if the list does not have any elements, return
  if(NULL==p)
  {
     return;
  }

  //the last element of the list is the only
  //element that has the "Next" field equal to NULL
  while(NULL!=p->Next)
  {
    pprev = p;
    p = p->Next;
  }
  
  //now "p" should give the last element
  //delete it
  gsl_vector_free(p->betaBayes);
  p->Next = NULL;
  delete p;

  //now the previous element in the list
  //becomes the last element
  pprev->Next = NULL;

  return;
}

//this function saves the regressions in the list with
//head "regressions" in a file with name "filename"
void SaveRegressions(char* filename,LPBayes regressions)
{
  int i;
  //open the output file
  FILE* out = fopen(filename,"w");
	
  if(NULL==out)
  {
    printf("Cannot open output file [%s]\n",filename);
    exit(1);
  }

  //this is the first regression
  LPBayes p = regressions->Next;
  while(NULL!=p)
  {
    //print the log marginal likelhood and the number of predictors
    fprintf(out,"Regression of Y on explanatory variable %d has",p->apredictor);
	fprintf(out,"log marginal likelihood %.5f (Laplace),  ",p->logmarglik);
	fprintf(out,"%.5f (Monte Carlo stable)",p->logmarglikMC);
	fprintf(out," with beta0 = %.5f ",gsl_vector_get(p->betaBayes,0));
	fprintf(out," with beta1 = %.5f ",gsl_vector_get(p->betaBayes,1));
    fprintf(out,"\n");

    //go to the next regression
    p = p->Next;
  }

  //close the output file
  fclose(out);

  return;
}
