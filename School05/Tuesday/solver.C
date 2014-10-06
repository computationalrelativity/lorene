#include "headcpp.h"
#include "math.h"

#include "tbl.h"
#include "matrice.h"
#include "leg.h"
#include "cheby.h"
		
		//*******************
		// First derivative
		//*******************
Tbl ope_der_cheb (const Tbl& so) {
	cout << "You have to work a bit : please implement me !" << endl ;
	abort() ;
	return so ;
}
		//*******************
		// Second derivative
		//*******************
Tbl ope_der_sec_cheb (const Tbl& so) {
	// Verification of size :
	assert (so.get_ndim() == 1) ;
	int size = so.get_dim(0) ;
	
	// The result is set to zero
	Tbl res (size) ;
	res.annule_hard() ;
	
	// the computation
	for (int n=0 ; n<size ; n++)
	     for (int p=n+2 ; p<size ; p++)
		if ((p+n)%2 == 0)
			res.set(n) += p*(p*p-n*n)*so(p) ;
	
	// Normalisation of the first coef :
	res.set(0) /= 2. ;
	
	return res ;
}


// Main code :

int main() {

     cout << "I am a lazy code... I dont do a thing ..." << endl ;
     return EXIT_SUCCESS ;
}
