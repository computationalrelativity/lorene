#include "headcpp.h"
#include "math.h"

#include "tbl.h"
#include "matrice.h"
#include "cheby.h"

double f_left (double) {
	return 1. ;
}

double f_right (double) {
	return 0. ;
}

double f_source (double x) {
	if (x<0)
	    return 1. ;
	else 
	   if (x>0) return 0. ;
	   else
	      return 0.5 ;
}

double f_solution (double x) {
	double ee = exp(1.) ;
	static double bmoins = -1./8./(1+ee*ee) - ee*ee/8./(1+ee*ee*ee*ee) ;
	static double bplus = ee*ee*ee*ee/8.*(ee*ee/(1+ee*ee*ee*ee)-1./(1+ee*ee)) ;
	
	if (x<0)
	    return 0.25 -(ee*ee/4.+bmoins*ee*ee*ee*ee)*exp(2*x)
	 		+ bmoins*exp(-2*x) ;
	else
	    return bplus*(exp(-2*x)-exp(2*x)/ee/ee/ee/ee) ;
}

		//*******************
		// First derivative
		//*******************
Tbl ope_der (const Tbl& so) {
	// Verification of size :
	assert (so.get_ndim() == 1) ;
	int size = so.get_dim(0) ;
	
	// The result is set to zero
	Tbl res (size) ;
	res.annule_hard() ;
	
	// the computation
	for (int n=0 ; n<size ; n++)
	     for (int p=n+1 ; p<size ; p++)
		if ((p+n)%2 == 1)
			res.set(n) += p*so(p)*2 ;
	
	// Normalisation of the first coef :
	res.set(0) /= 2. ;
	
	return res ;
}
		//*******************
		// Second derivative
		//*******************
Tbl ope_der_sec (const Tbl& so) {
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

		//**************************************
		//    The operator (multi domain case)
		//**************************************

Matrice multi_ope (int n) {

	// The result :
	Matrice res(n,n) ;
	res.set_etat_qcq() ;
	
	// Work arrays :
	Tbl so (n) ;
	Tbl dd_so (n) ;
	
	// Column by column :
	for (int col=0 ; col<n ; col++) {
		so.annule_hard() ;
		so.set(col) = 1 ;
		
		// The derivative
		dd_so = ope_der_sec(so) ;
		
		// Put in the matrix :
		for (int line=0 ; line<n ; line++)
		    res.set(line, col) = -4*dd_so(line) + 4*so(line) ;
	}

	return res ;
}

int main () {

	int nr ;
	cout << "Please enter the number of coefficients :" << endl ;
	cin >> nr ;
	
	cout << "The operator matrix is " << endl ;
	Matrice operator_mat (multi_ope(nr)) ;
	cout << operator_mat << endl ;
	
	// Resolution with a tau method :
	Matrice systeme(nr*2, nr*2) ;
	systeme.annule_hard() ;
	
	// Boundary conditions are inforced by additional equations:	
	// Left boundary condition :
	for (int i=0 ; i<nr ; i++)
	    systeme.set(0, i) = (i%2==0) ? 1 : -1 ;
	// Equation in the first domain :
	for (int i=0 ; i<nr-2 ; i++)
	    for (int j=0 ; j<nr ; j++)
	        systeme.set(i+1,j) = operator_mat(i,j) ;
	// Continuity of the solution :
	for (int i=0 ; i<nr ; i++)
	    systeme.set(nr-1, i) = 1 ;
	for (int i=0 ; i<nr ; i++)
	    systeme.set(nr-1, i+nr) = (i%2==0) ? -1 : 1 ;
	    
	// Continuity of the first derivative :
	for (int i=0 ; i<nr ; i++)
	    systeme.set(nr, i) = i*i ;
	for (int i=0 ; i<nr ; i++)
	    systeme.set(nr, i+nr) = (i%2==0) ? i*i : -i*i ;
	
	// Equation in the second domain :
	for (int i=0 ; i<nr-2 ; i++)
	    for (int j=0 ; j<nr ; j++)
	        systeme.set(i+1+nr,j+nr) = operator_mat(i,j) ;
	    
	// Right boundary condition :
	for (int i=0 ; i<nr ; i++)
	    systeme.set(2*nr-1, i+nr) = 1 ;
	    
	systeme.set_lu() ;
	
	cout << "Multi-domain systeme of equations :" << endl ;
	cout << systeme << endl ;
	    
	// Coeficients of the source :
	// left hand side :
	Tbl conf (nr) ;
	conf.set_etat_qcq() ;
	Tbl colocation(coloc_cheb(nr)) ;
	for (int i=0 ; i<nr ; i++)
	    conf.set(i) = f_left(0.5*(colocation(i)-1)) ;
	Tbl coefs_left (coef_cheb(conf)) ;
	for (int i=0 ; i<nr ; i++)
	    conf.set(i) = f_right(0.5*(colocation(i)+1)) ;
	Tbl coefs_right (coef_cheb(conf)) ;
	
	
	cout << "Coefficients of the source, left and right" << endl ;
	cout << coefs_left << endl ;
	cout << coefs_right << endl ;
	
	// Second member of the system :
	Tbl sec_member (2*nr) ;
	sec_member.set_etat_qcq() ;
	// Leftoundary conditions :
	sec_member.set(0) = 0 ;
	// Source on the left:
	for (int i=1 ; i<nr-1 ; i++)
	    sec_member.set(i) = coefs_left(i-1) ;
	// continuity of the function :
	sec_member.set(nr-1) = 0 ;
	// continuity of the first derivative :
	sec_member.set(nr) = 0 ;
	// Source on the right:
	for (int i=1 ; i<nr-1 ; i++)
	    sec_member.set(i+nr) = coefs_right(i-1) ;
	// Right boundary condition :
	sec_member.set(2*nr-1) = 0 ;
	
	cout << "Second member for the multi-domain system" << endl ;
	cout << sec_member << endl ;
	
	// The system is inverted :
	Tbl coef_sol (systeme.inverse(sec_member)) ;
	    
	cout << "Coefficients of the solution : " << endl ;
	cout << coef_sol << endl ;
	
	// Output in a file for plotting purposes :
	char name_out[30] ;
	sprintf(name_out, "plot_multi_tau_%i.dat", nr) ;
	ofstream fiche (name_out) ;
	
	int resolution = 200 ;
	double x=-1 ;
	double step = 2./resolution ;
	
	// We will also compute the maximum difference between the solution and its numerical value :
	double error_max = 0 ;
	
	double xi ;
	int offset ;
	for (int i=0 ;  i<resolution+1 ; i++)  {
	
	    // in which domain ?
	    if (x<0) {
	        offset=0 ;
		xi = 2*x+1 ;
		}
	    else {
	        offset=nr ;
		xi = 2*x-1 ;
	    }
	    
	    // One computes the solution at the current point
	    double val_func = 0 ;
	    for (int j=0 ; j<nr ; j++)
	         val_func += coef_sol(j+offset)*cheby(j, xi) ;
	   fiche << x << " " << val_func << " " << f_solution(x) << endl ;
	   double error = fabs(val_func-f_solution(x)) ;
	   if (error > error_max)
	        error_max = error ;
	   
	   x += step ;
	}
	
	cout << "Error max : " << endl ;
	cout << error_max << endl ;
	
	return 0 ;
}
