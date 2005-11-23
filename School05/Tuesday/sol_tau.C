#include "headcpp.h"
#include "math.h"

#include "tbl.h"
#include "matrice.h"
#include "cheby.h"

double f_source (double x) {
	static double constante = -4*exp(1.)/(1+exp(1.)*exp(1.)) ;
	return constante + exp(x) ;
}

double f_solution (double x) {
	static double constante = -4*exp(1.)/(1+exp(1.)*exp(1.)) ;
	return exp(x) - sinh(1.)/sinh(2.)*exp(2*x) + constante/4 ;
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

		//*******************
		//    The operator
		//*******************

Matrice matrix_ope (int n) {

	// The result :
	Matrice res(n,n) ;
	res.set_etat_qcq() ;
	
	// Work arrays :
	Tbl so (n) ;
	Tbl d_so (n) ;
	Tbl dd_so (n) ;
	
	// Column by column :
	for (int col=0 ; col<n ; col++) {
		so.annule_hard() ;
		so.set(col) = 1 ;
		
		// The derivatives
		d_so = ope_der(so) ;
		dd_so = ope_der_sec(so) ;
		
		// Put in the matrix :
		for (int line=0 ; line<n ; line++)
		    res.set(line, col) = dd_so(line) -4*d_so(line) + 4*so(line) ;
	}

	return res ;
}
int main () {

	int nr ;
	cout << "Please enter the number of coefficients :" << endl ;
	cin >> nr ;
	
	cout << "The operator matrix is " << endl ;
	Matrice operator_mat (matrix_ope(nr)) ;
	cout << operator_mat << endl ;
	
	// Resolution with a tau method :
	Matrice systeme(nr, nr) ;
	systeme.set_etat_qcq() ;
	// Boundary conditions are inforced by additional equations:	
	// Left boundary condition :
	for (int i=0 ; i<nr ; i++)
	    systeme.set(0, i) = (i%2==0) ? 1 : -1 ;
	// Right boundary condition :
	for (int i=0 ; i<nr ; i++)
	    systeme.set(1, i) = 1 ;
	// The rest of the equations :
	for (int i=2 ; i<nr ; i++)
	    for (int j=0 ; j<nr ; j++)
	        systeme.set(i,j) = operator_mat(i-2, j) ;

	systeme.set_lu() ;
		
	cout << "The Tau matrix is" << endl ;
	cout << systeme << endl ;
	
	// Coeficients of the source :
	Tbl conf (nr) ;
	conf.set_etat_qcq() ;
	Tbl colocation(coloc_cheb(nr)) ;
	for (int i=0 ; i<nr ; i++)
	    conf.set(i) = f_source(colocation(i)) ;
	Tbl coefs (coef_cheb(conf)) ;
	
	cout << "Coefficients of the source" << endl ;
	cout << coefs << endl ;
	
	// Second member of the system :
	Tbl sec_member (nr) ;
	sec_member.set_etat_qcq() ;
	// Boundary conditions :
	sec_member.set(0) = 0 ;
	sec_member.set(1) = 0 ;
	// Coefficients of the source :
	for (int i=2 ; i<nr ; i++)
	    sec_member.set(i) = coefs(i-2) ;
	
	cout << "Second member for the Tau system" << endl ;
	cout << sec_member << endl ;
	   
	// The system is inverted :
	Tbl coef_sol (systeme.inverse(sec_member)) ;
	    
	cout << "Coefficients of the solution : " << endl ;
	cout << coef_sol << endl ;
	
	// Output in a file for plotting purposes :
	char name_out[30] ;
	sprintf(name_out, "plot_tau_%i.dat", nr) ;
	ofstream fiche (name_out) ;
	
	int resolution = 200 ;
	double x=-1 ;
	double step = 2./resolution ;
	
	// We will also compute the maximum difference between the solution and its numerical value :
	double error_max = 0 ;
	
	for (int i=0 ;  i<resolution+1 ; i++)  {
	    // One computes the solution at the current point
	    double val_func = 0 ;
	    for (int j=0 ; j<nr ; j++)
	         val_func += coef_sol(j)*cheby(j, x) ;
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
