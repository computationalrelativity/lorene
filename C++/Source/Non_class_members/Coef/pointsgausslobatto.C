#include <math.h>
#include "tbl.h"
#include <cassert>

namespace Lorene {
double* jacobi(int , double) ;

double* pointsgausslobatto(int n) {

  int nmax = 10 ;
  int ndiv = (n+1)*(n+1)+5 ;
  double pas = double(2)/double(ndiv-1) ;
  double eps = 2.5E-12 ;

  Tbl subdiv(ndiv) ; 
  subdiv.set_etat_qcq();
  Tbl cs(n-1,2) ;
  cs.set_etat_qcq() ;
  double* xx = new double[nmax] ;
  double* yy ;
  double* zz ;
  int i,k,l ;

  for (i = 0 ; i < ndiv ; i++) {
    subdiv.set(i) = double(-1) + pas * double(i) ;
  }

  double a = (2*double(n)+3)/double((n+1)*(n+1)) ;
  double b = - 1 - a ;

  int j=0;

  for (i = 1 ; i < ndiv-3 ; i++) {
    yy = jacobi(n+1 , subdiv(i)) ;
    zz = jacobi(n+1 , subdiv(i+1)) ;
    double omega1 = yy[n+1] + a * yy[n] + b * yy[n-1] ;
    double omega2 = zz[n+1] + a * zz[n] + b * zz[n-1] ;
    if (omega1*omega2 <= 0) {
      cs.set(j,0) = subdiv(i) ;
      cs.set(j,1) = subdiv(i+1) ;
      j++;
    }
    delete [] yy ;
    delete [] zz ;
}


  
  double* pointsgl = new double[j+2];
  assert(j==n-1) ;
  pointsgl[0] = -1;

  for (l = 0 ; l < j ; l++) {
    xx[0] = cs(l,0) ; 
    xx[1] = cs(l,1) ;
    for (k = 2 ; k < nmax ; k++) {
       yy = jacobi(n+1 , xx[k-2]) ;
       zz = jacobi(n+1 , xx[k-1]) ;
       double omega1 = yy[n+1] + a * yy[n] + b * yy[n-1] ;
       double omega2 = zz[n+1] + a * zz[n] + b * zz[n-1] ;
       if (( fabs(xx[k-2]-xx[k-1])>=eps )&&(fabs(omega2)>=eps )) {
	 xx[k]=(xx[k-2]*omega2 -xx[k-1]*omega1)/(omega2-omega1) ;
       }
       else {
	 xx[k]=xx[k-1];
       }
       delete [] yy ;
       delete [] zz ;
    }
    pointsgl[l+1] = xx[nmax-1] ;
  }
  delete [] xx ;
  
  pointsgl[n] = 1 ;
	
  return pointsgl ;
}
}
