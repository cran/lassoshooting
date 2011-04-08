#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <math.h>

#ifndef CBLAS
#  include <R_ext/BLAS.h>
#  define daxpy F77_CALL(daxpy)
#else
#  include <cblas.h>
#endif

#include "ccd.h"

inline double softthresh(double x,double t) {
  double v = fabs(x) - t;
  double ret = 0.;
  if (v  > 0.) {
    if (x >= 0.) ret = v;
    else ret = -v;
  } 
//  printf("softthresh: %f\n", ret);
	return ret;
}
int ccd_common(param_t* params) {
  int its = 0;
  double delta = 0.0;
  double deltabeta = 0.0;
  double betajstar = 0.0;
  double betajstar_old = 0.0;
  //bool penalizethis = true;
  //int state = 1;
  int p = params->p;
  double *s = params->Xty;
  double factor2 = params->factor2;

  // inf-norm of X'y
  double infnorm = 0.;
  for (int i=0; i < p; ++i) {
    double this = fabs(s[i] / params->factor2);
    if (this > infnorm)
      infnorm = this;
  }
  if (params->trace > 0) printf("lambda: %f\n",params->lambda);
  if (params->trace > 0) printf("infnorm: %f\n",infnorm);
  if (params->lambda > infnorm) {
    if (params->trace > 0) printf("returning because lambda > infnorm\n");
    return 1; // XXX: quit if lambda > ||X'y||_inf
  }


  if (params->trace >= 3)
    for (int i=0; i < p; ++i) {
      printf("penalize beta_%d with %.2f\n",i,params->w[i]);
    }
  //

  do {
    delta = 0.0;
    for(int j=0;j < p; ++j) {
      //gsl_vector_view XtXj = gsl_matrix_column(XtX_, (size_t)j);
      //double XtXjj = gsl_matrix_get(XtX_,j,j);
      //double XtXjj = gsl_vector_get(&XtXj.vector, j);
      double XtXjj = params->XtX[j+j*p];
//			printf("XtXjj = %f  \n",XtXjj);

      if (XtXjj == 0. || params->forcezero == j+1) continue;
      //if (active && its % 10 != 0 && beta_[j] == 0.)  continue;
 
      betajstar_old = betajstar;
      betajstar = s[j] + (factor2*XtXjj) * params->beta[j];
      if (params->trace >= 2) {
	FILE*D = fopen("ccd.debug","a");
	fprintf(D,"forcezero=%d, its=%d, betajstar-pre: %f\n", params->forcezero,its,betajstar);
	fclose(D);
      }

      if (isinf(betajstar) || isnan(betajstar)) {
	fprintf(stderr,"******************************************\n"
	     __FILE__ ": BUG OR PATHOLOGICAL DATA\n\n");
	fprintf(stderr, "Please mail me the data that can reproduce this error <Tobias.Abenius@Chalmers.SE>\n");
	fprintf(stderr, "betajstar prev = %f  \n",betajstar_old);
	fprintf(stderr, "deltabeta prev = %f  \n",deltabeta);
	fprintf(stderr, "s_%d = %f  \n",j,s[j]);
	fprintf(stderr, "betajstar_%d = %f  \n",j,betajstar);
	fprintf(stderr, "beta_%d = %f  \n",j,params->beta[j]);
	fprintf(stderr, "XtXjj = %f  \n",XtXjj);
	fprintf(stderr, "\nGiving up...\n");
	fprintf(stderr,"******************************************\n");
	return 0;
      }
//      if (fabs(params->w[j]*params->lambda) < 1e-40) {
//        betajstar /= factor2*XtXjj;
//      } else {
        betajstar = softthresh(betajstar, params->w[j]*params->lambda) / (factor2*XtXjj);
      if (params->trace >= 2) {
	FILE*D = fopen("ccd.debug","a");
	fprintf(D,"forcezero=%d, its=%d, betajstar-post: %f, lam: %f, div: %f\n", params->forcezero,its,betajstar,params->w[j]*params->lambda,factor2*XtXjj);
	fclose(D);
      }

//      }
			//printf("beta_%d = %.4f\n",j,betajstar);
      deltabeta = betajstar - params->beta[j];
      if (params->trace >= 2) {
	FILE*D = fopen("ccd.debug","a");
	fprintf(D,"forcezero=%d, its=%d, deltabeta: %f\n", params->forcezero,its,deltabeta);
	fclose(D);
      }

      params->beta[j] = betajstar;
      delta = max(delta, fabs(deltabeta));
      /* s <- s - 2*deltabeta XtX(:,j);
         axpy :: y <- ax + y*/
#ifndef CBLAS
#warning Using R fortran BLAS calls
      const double factor = -deltabeta*factor2;
      int one = 1; // memory layout
      daxpy(&p, &factor, &params->XtX[j], &p, s, &one);
#else
#warning Using cblas BLAS calls
      cblas_daxpy(p, -deltabeta*factor2,&params->XtX[j],p, s,1);
#endif
      if (params->trace >= 2) 
	for (int di=0; di < p; ++di)  {
	  FILE*D = fopen("ccd.debug","a");
	  fprintf(D,"forcezero=%d, its=%d, s%d: %f\n", params->forcezero,its,di,s[di]);
	  fclose(D);
	}

      //gsl_blas_daxpy(-deltabeta*factor2,&XtXj.vector, &s.vector);
    }
    if (params->trace >= 2) {
      FILE*D = fopen("ccd.debug","a");
      fprintf(D,"its = %d \tdelta %f  \n",its,delta);
      fclose(D);
    }
  } while (++its < params->maxits && delta > params->tol);
  if (params->trace) printf("ccd ran for %d iterations, delta: %g\n",its,delta);
  params->its = its;
  params->delta = delta;
  for(int i=0; i < params->p; ++i) {
    params->beta[i] *= factor2;
  }
  return 1;
}

