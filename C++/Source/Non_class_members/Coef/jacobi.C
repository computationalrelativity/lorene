namespace Lorene {

double* jacobi(int n, double x) {

  int i ;
  double* J = new double[n+1] ;

  if (n==0) {

    J[0] = double(1) ;
  }
  else {

    J[0] = double(1) ;
    J[1] = double(2) * x - 1 ;
    for ( i = 2 ; i < n+1 ; i++) {
      double l = double(i) ;
      J[i] = ((2*l + 1)*(l*(l+1)*x - 1)*J[i-1] - (l-1)*(l+1)*(l+1)*J[i-2])/(l*l*(l+2)) ;

    }
  }
  return J ;
}
}
