#include "headcpp.h"
#include "math.h"

#include "tbl.h"
#include "matrice.h"
#include "cheby.h"

// Chebyshev polynomial
double cheby (int n, double x) {

	if (n==0)
	    return 1 ;
	else if (n==1)
	    return x ;
	     else {
	         double tm2 = 1 ;
		 double tm1 = x ;
		 double tm ;
		 for (int i=2 ; i<=n ; i++) {
		     tm = 2*x*tm1 - tm2 ;
		     tm2 = tm1 ;
		     tm1 = tm ;
		 }
		 return tm ;
	}
}

// Norme
double norme_cheb (int i) {
    if (i==0)
        return M_PI ;
    else
        return M_PI/2 ;
}

// Colocation point
Tbl coloc_cheb (int n) {
     Tbl res (n) ;
     res.set_etat_qcq() ;
     for (int i=0 ; i<n ; i++)
     res.set(i) = -cos(M_PI*i/(n-1)) ;
     return res ;
}
 
// Weight functions
Tbl weight_cheb (int n) {
    Tbl res (n) ;
    res.set_etat_qcq() ;
    for (int i=0 ; i<n ; i++)
    	res.set(i) = ((i==0) || (i==n-1)) ? M_PI/2./n : M_PI/n ;
    return res ;
}

// Normalisation factor :
Tbl gamma_cheb (int n) {
	Tbl res (n) ;
	res.annule_hard() ;
	for (int i=0 ; i<n ; i++)
 	   	res.set(i) = M_PI/2*(n-1)/n ;
	res.set(0) *= 2 ;
	res.set(n-1) *= 2 ;
	return res ;
}

// Coefficients :
Tbl coef_cheb (const Tbl& so) {
	
	// Verification of size :
	assert (so.get_ndim() == 1) ;
	int n = so.get_dim(0) ;
	
	// Auxiliary quantities :
	Tbl colocation (coloc_cheb(n)) ;
	Tbl poids (weight_cheb(n)) ;
	Tbl norme (gamma_cheb(n)) ;
	
	// The result is set to zero
	Tbl res (n) ;
	res.annule_hard() ;

	for (int i=0 ; i<n ; i++)
	    for (int j=0 ; j<n ; j++) 
	     	res.set(i) += so(j)*cheby(i, colocation(j))*poids(j)/norme(i) ;
	return res ;
}
