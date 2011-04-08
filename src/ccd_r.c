#define HAVE_INLINE
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <math.h>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h>	
#include <R_ext/BLAS.h>

#include "ccd.h"

SEXP ccd(SEXP args) {
  param_t params;
  double *X,*y;
  double *givenXtX=NULL,*givenXty=NULL;
  int Xm, Xn;
  int Ym, Yn;
  int XtXp,Xtym;
  Xtym = XtXp = Xm = Xn = Ym = Yn = -1;
  int wm = -1,wn = -1;

  double factor2=2.;

  //default parameters
  params.m = params.p = -1;
  params.tol = 1e-6;
  params.forcezero = -1;
  params.maxits = 10000;
  params.trace = 0;
  int *nopenalize = NULL;
  params.w = NULL;
  params.factor2 = factor2;

  // argument handling
  args = CDR(args); 
  for(int i = 0; args != R_NilValue; i++, args = CDR(args)) {
    const char *name = CHAR(PRINTNAME(TAG(args)));
    if (params.trace > 1) {
      Rprintf(__FILE__ ": parsing parameter %s\n",name);
    }
    if (CAR(args) == R_NilValue) {
      //if (params.trace) Rprintf(__FILE__ ": parameter %s was null\n",name);
      continue;
    }
    if (strcasecmp(name, "x")==0) { 
      X = REAL(CAR(args)); 
      SEXP d = getAttrib(CAR(args), R_DimSymbol); // dimensions
      Xm = INTEGER(d)[0];
      Xn = INTEGER(d)[1];
    }
    else if (strcasecmp(name, "y")==0) { 
      y = REAL(CAR(args)); 
      SEXP d = getAttrib(CAR(args), R_DimSymbol); // dimensions
      if(!isNull(d)) {
        Ym = INTEGER(d)[0];
        Yn = INTEGER(d)[1];
      } else {
        Ym = length(CAR(args));
        Yn = 1;
      }
    }
    else if (strcasecmp(name, "lambda")==0)  { 
      if (length(CAR(args)) != 1) { error("length of lambda should be 1!\n"); }
      params.lambda = REAL(CAR(args))[0]; 
    }
    else if (strcasecmp(name, "forcezero")==0) { 
      if (length(CAR(args)) != 1) { error("length of forcezero should be 1!\n"); }
      SEXP newforcezero;
      //double* nopen = REAL(CAR(args)); 
      PROTECT(newforcezero = AS_INTEGER(CAR(args)));
      params.forcezero = INTEGER(newforcezero)[0]; 
      UNPROTECT(1);
    } 
    else if (strcasecmp(name, "thr")==0) { 
      if (length(CAR(args)) != 1) { error("length of thr should be 1!\n"); }
      params.tol = REAL(CAR(args))[0]; 
    } 
    else if (strcasecmp(name, "maxit")==0) { 
      if (length(CAR(args)) != 1) { error("length of maxit should be 1!\n"); }
      params.maxits = REAL(CAR(args))[0]; 
    } 
    else if (strcasecmp(name, "penaltyweight")==0) { 
      params.w = REAL(CAR(args)); 
      SEXP d = getAttrib(CAR(args), R_DimSymbol); // dimensions
      if(!isNull(d)) {
        wm = INTEGER(d)[0];
        wn = INTEGER(d)[1];
      } else {
        wm = length(CAR(args));
        wn = 1;
      }
      if (wn != 1) error("penaltyweight should be p x 1");
    }
    else if (strcasecmp(name, "nopenalize")==0) { 
      SEXP newnopen;
      //double* nopen = REAL(CAR(args)); 
      PROTECT(newnopen = AS_NUMERIC(CAR(args)));
      double *nopen = REAL(newnopen);
      UNPROTECT(1);
      int N = length(CAR(args));
      nopenalize = (int*)R_alloc(N+1,sizeof(int));
      for(int i=0; i < N; ++i) {
        nopenalize[i] = nopen[i];
        if (nopenalize[i] < 0) 
          error("Can not penalize variables with negative index! The range is 0 to p-1.\n"); 
      }
      nopenalize[N] = -1;
    } 
    else if (strcasecmp(name, "trace")==0) { 
      if (length(CAR(args)) != 1) { error("length of trace should be 1!\n"); }
      params.trace = REAL(CAR(args))[0]; 
      if (params.trace > 0) Rprintf("Tracing on!\n");
    } 
    else if (strcasecmp(name,"XtX")==0) {
      givenXtX = REAL(CAR(args)); 
      SEXP d = getAttrib(CAR(args), R_DimSymbol); // dimensions
      int m = INTEGER(d)[0];
      int n = INTEGER(d)[1];
      if (m!=n) error("X'X should be square, it is %dx%d!",m,n);
      XtXp = m;
    }
    else if (strcasecmp(name,"Xty")==0) {
      givenXty = REAL(CAR(args)); 
      SEXP d = getAttrib(CAR(args), R_DimSymbol); // dimensions
      int m = INTEGER(d)[0];
      int n = INTEGER(d)[1];
      if (n!=1) error("X'y should be Mx1, it is %dx%d!",m,n);
      Xtym = m;
    }
    else {
      error("Unknown parameter '%s'!\n",name); 
      return(R_NilValue);
    }
  }
  if (params.trace >=2) {
    Rprintf("X m: %d n: %d\n",Xm,Xn); 
    Rprintf("Y m: %d n: %d\n",Ym,Yn); 
  }
  int p = (Xn == -1)? XtXp : Xn; // no variables
  int m = (Xm == -1)? Xtym : Xm; // no equations
  if (wm != -1 && wm != p) error("penaltyweight should be p x 1\n");
  params.p = p;
  params.m = m;
  if (!params.w) {
    params.w = (double*)R_alloc(p,sizeof(double));
    for(int i=0; i < p; ++i) params.w[i] = 1.0f;
  }
  for(int i=0; nopenalize && nopenalize[i] >= 0; ++i) {
    params.w[nopenalize[i]] = 0.0f;
  }

  int ret;
  params.XtX = (double*)R_alloc(p*p,sizeof(double));

  double zero=0.;
  double one=1;
  if (params.trace>=2)  printf("havextx: %d\n",givenXtX!=NULL);
  if (!givenXtX && Xm != -1 && Xn != -1) { 
    if (params.trace>=2) printf("calculating X'X\n",factor2);
    // X is MxN   X'X  NxM*MxN -> NxN
    F77_CALL(dgemm)("T","N",&Xn,&Xn,&Xm, &one,X,&Xm, X,&Xm, &zero,params.XtX,&Xn);  // 2X'X -> XtX
  } 
  else if (givenXtX) {
    if (params.trace>=2) printf("givenXtX\n");
    for(int i=0; i < p*p; ++i)
      params.XtX[i] = givenXtX[i];
  }
  else {
    error("Need either X and y or XtX and Xty");
  }
  if (params.trace>=2) printf("scaling X'X with %f, p: %d\n",factor2,p);
  for(int i=0; i < p*p; ++i)
    params.XtX[i] = factor2*params.XtX[i];
  if (params.trace >= 2) {
    FILE*D = fopen("ccd.debug","a");
    for (int i=0; i < p; ++i)
      for (int j=0; j < p; ++j)
	fprintf(D,"X'X[%d,%d]: %f\n", i,j,params.XtX[i*p+j]);
    fclose(D);
  }


  params.Xty = (double*)R_alloc(p,sizeof(double)); // NxM*Mx1 -> Nx1
  if (!givenXty && Xm != -1 && Ym != -1) {
    F77_CALL(dgemv)("T", &Xm,&Xn, &factor2,X,&Xm, y,&Yn, &zero,params.Xty,&Yn); // Xty <- 2X'y
  } else if (givenXty) {
    for(int i=0; i < Xtym; ++i)
      params.Xty[i] = factor2*givenXty[i];
  } else {
    error("Need either X and Y or XtX and Xty");
  }
  if (params.trace > 2) printf("X'y_(1,1)=%.4f\n",params.Xty[0]/2.);
  if (params.trace > 2) printf("X'y_(2,1)=%.4f\n",params.Xty[1]/2.);
  if (params.trace > 2) printf("X'y_(3,1)=%.4f\n",params.Xty[2]/2.);

  SEXP beta_;
  PROTECT(beta_ = allocVector(REALSXP,p));
  params.beta = REAL(beta_);
  for(int i=0; i < params.p; ++i) {
    params.beta[i] = 0.;
  }
  ccd_common(&params);

#define NC 3
  char* names[NC] = {"coefficients","iterations","delta"};
  SEXP list_names, list, itsR, deltaR;
  PROTECT(itsR = NEW_INTEGER(1));
  int* itsRp = INTEGER_POINTER(itsR);
  *itsRp = params.its;
  PROTECT(deltaR = NEW_NUMERIC(1));
  double* deltaRp = NUMERIC_POINTER(deltaR);
  *deltaRp = params.delta;

  PROTECT(list_names = allocVector(STRSXP, NC));
  for(int i = 0; i < NC; i++)
    SET_STRING_ELT(list_names, i,  mkChar(names[i]));
 
  PROTECT(list = allocVector(VECSXP, NC)); // Creating a list with NC vector elements
  SET_VECTOR_ELT(list, 0, beta_);         // attaching beta vector to list
  SET_VECTOR_ELT(list, 1, itsR);      // attaching its vector to list
  SET_VECTOR_ELT(list, 2, deltaR);      // attaching its vector to list
  setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names
  UNPROTECT(5);
  return list;
  /*
  p = 3;
  PROTECT(beta_ = allocVector(REALSXP, p));
  double *beta = REAL(beta_);
  beta[0] = 1;
  beta[1] = 3;
  beta[2] = 5;
  UNPROTECT(1);
  return beta_;
  */
} 

