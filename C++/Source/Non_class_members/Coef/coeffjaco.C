#include "tbl.h"

Tbl jacobipointsgl(int) ;

double* coeffjaco(int n , double* ff) {

  Tbl jj = jacobipointsgl(n) ;
  double som ;
  int i,k ;
  double* aa = new double[n+1] ;

  for (k = 0 ; k < n ; k++) {
    som = 3*ff[0]*jj(k,0)/(jj(n,0)*jj(n,0)) ;
    for (i = 1 ; i < n+1 ; i++) {
      som += ff[i]*jj(k,i)/(jj(n,i)*jj(n,i)) ;
    }
    aa[k] = (2*k+3)/double(n*(n+3))*som ;
  }
  som = 3*ff[0]/jj(n,0) ;
    for (i = 1 ; i < n+1 ; i++) {
    som += ff[i]/jj(n,i) ;
    }
    aa[n]=som/double(n+3) ;
  return aa ;
}
