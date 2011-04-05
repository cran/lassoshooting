#ifndef _CCD_H_
#  define _CCD_H_ 1

#  define max(a,b) (((a)>(b))?(a):(b))

typedef struct params {
  double *XtX;
  double *Xty;
  double lambda;
  double *beta;
  int m;
  int p;

  int forcezero;
  int maxits;
  int its;
  double delta;
  double tol;
  int trace;
  double *w;
  double factor2;
} param_t;

double softthresh(double x,double t);
int ccd_common(param_t* params);

#endif
