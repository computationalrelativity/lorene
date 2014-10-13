#include "tbl.h"

namespace Lorene {
double* jacobi(int , double) ;
double* pointsgausslobatto(int) ;

Tbl jacobipointsgl(int n) {

 double* pointsgl = pointsgausslobatto(n) ;

 Tbl jj(n+1,n+1) ;
 jj.set_etat_qcq() ;

 int i,k ;

 for (i = 0 ; i < n+1 ; i++ ) {
    double* yy = jacobi(n,pointsgl[i]) ;
   for (k = 0 ; k < n+1 ; k++ ) {
     jj.set(k,i) = yy[k] ;
   }
   delete [] yy ;
 }
 delete [] pointsgl ;
 return jj ;
}
}
