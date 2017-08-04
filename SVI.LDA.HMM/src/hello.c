#include <R.h>
#include <Rmath.h>
#include<Rinternals.h>
SEXP vecSum(SEXP Rvec)
{
  int i,n;
  double *vec, value = 0;
  vec = REAL(Rvec);
  n = length(Rvec);
  for(i=0; i< n;i++)
  {
    value += vec[i];
  }
  printf("The total value is: %4.6f \n", value);
  return R_NilValue;
}
SEXP ab(SEXP Ra, SEXP Rb){
  int i, a, b;
  SEXP Rval;
  Ra = coerceVector(Ra, INTSXP);
  Rb = coerceVector(Rb, INTSXP);
  a = INTEGER(Ra)[0];
  b = INTEGER(Rb)[0];
  PROTECT(Rval = allocVector(INTSXP, b - a + 1));
  for (i = a; i <= b; i++)
    INTEGER(Rval)[i - a] = i;
  UNPROTECT(1);
  return Rval;
}
