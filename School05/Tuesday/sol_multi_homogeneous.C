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
	
	// The non degenerare operator :
	Matrice nondege (nr-2, nr-2) ;
	nondege.set_etat_qcq() ;
	for (int lig=0 ; lig<nr-2 ; lig++)
	    for (int col=0 ; col<nr-2 ; col++)
	         nondege.set(lig, col) = operator_mat(lig, col+2) ;
	nondege.set_lu() ;
	cout << "The inverted matrix is :" << endl ;
	cout << nondege << endl ;
	
	// Computation of the particular solution on the left :
	Tbl colocation(coloc_cheb(nr)) ;
	Tbl so_left (nr) ;
	so_left.set_etat_qcq() ;
	for (int i=0 ; i<nr ; i++)
	     so_left.set(i) = f_left (0.5*(colocation(i)-1)) ;
	Tbl coef_so_left (coef_cheb(so_left)) ;
	
	cout << "Coefficients of the source on the left : " << endl ;	
	cout << coef_so_left << endl ;
	
	Tbl auxi_so_left (nr-2) ;
	auxi_so_left.set_etat_qcq() ;
	for (int i=0 ; i<nr-2 ; i++)
		auxi_so_left.set(i) = coef_so_left(i) ;
	Tbl inv_left (nondege.inverse(auxi_so_left)) ;
	Tbl sp_left (nr) ;
	sp_left.set_etat_qcq() ;
	sp_left.set(0) = 0 ;
	sp_left.set(1) = 0 ;
	for (int i=0 ; i<nr-2 ; i++)
	    sp_left.set(i+2) = inv_left(i) ;
	 
	cout << "Particular solution on the left :" << endl ;
	cout << sp_left << endl ;
	
	// Computation of the particular solution on the right :
	Tbl so_right (nr) ;
	so_right.set_etat_qcq() ;
	for (int i=0 ; i<nr ; i++)
	     so_right.set(i) = f_right (0.5*(colocation(i)+1)) ;
	Tbl coef_so_right (coef_cheb(so_right)) ;
	
	cout << "Coefficients of the source on the right : " << endl ;	
	cout << coef_so_right << endl ;
	
	Tbl auxi_so_right (nr-2) ;
	auxi_so_right.set_etat_qcq() ;
	for (int i=0 ; i<nr-2 ; i++)
		auxi_so_right.set(i) = coef_so_right(i) ;
	Tbl inv_right (nondege.inverse(auxi_so_right)) ;
	Tbl sp_right (nr) ;
	sp_right.set_etat_qcq() ;
	sp_right.set(0) = 0 ;
	sp_right.set(1) = 0 ;
	for (int i=0 ; i<nr-2 ; i++)
	    sp_right.set(i+2) = inv_right(i) ;
	 
	cout << "Particular solution on the right :" << endl ;
	cout << sp_right << endl ;
	
	// Computation of the homogeneous solutions :
	Tbl val_sh_plus_left (nr) ;
	val_sh_plus_left.set_etat_qcq() ;
	for (int i=0 ; i<nr ; i++)
	     val_sh_plus_left.set(i) = exp(colocation(i)-1) ;
	 Tbl coef_sh_plus_left (coef_cheb(val_sh_plus_left)) ;
	 
	 Tbl val_sh_moins_left (nr) ;
	 val_sh_moins_left.set_etat_qcq() ;
	for (int i=0 ; i<nr ; i++)
	     val_sh_moins_left.set(i) = exp(-colocation(i)+1) ;
	 Tbl coef_sh_moins_left (coef_cheb(val_sh_moins_left)) ;
	 
	Tbl val_sh_plus_right (nr) ;
	val_sh_plus_right.set_etat_qcq() ;
	for (int i=0 ; i<nr ; i++)
	     val_sh_plus_right.set(i) = exp(colocation(i)+1) ;
	 Tbl coef_sh_plus_right (coef_cheb(val_sh_plus_right)) ;
	 
	 Tbl val_sh_moins_right (nr) ;
	 val_sh_moins_right.set_etat_qcq() ;
	 for (int i=0 ; i<nr ; i++)
	     val_sh_moins_right.set(i) = exp(-colocation(i)-1) ;
	 Tbl coef_sh_moins_right (coef_cheb(val_sh_moins_right)) ;
	 
	 cout << "Various homogeneous solutions : " << endl ;
	 cout << coef_sh_plus_left << endl ;
	 cout << coef_sh_moins_left << endl ;
	 cout << coef_sh_plus_right << endl ;
	 cout << coef_sh_moins_right << endl ;
	 
	 // Construction of the matching matrix :
	 Matrice matching (4,4) ;
	 matching.annule_hard() ;
	 // Boundary condition on the left :
	 matching.set(0,0) = exp(2.) ;
	 matching.set(0,1) = exp(-2.) ;
	 // matching of the solution accross the boundary :
	 matching.set(1,0) = 1 ;
	 matching.set(1,1) = 1 ;
	 matching.set(1,2) = -1 ;
	 matching.set(1,3) = -1 ; 
	// matching of the first derivative :
	matching.set(2,0) = 1 ;
	matching.set(2,1) = -1 ;
	matching.set(2,2) = 1 ;
	matching.set(2,3) = -1 ;
	// boundary on the right :
	matching.set(3,2) = exp(2.) ;
	matching.set(3,3) = exp(-2.) ;
	
	matching.set_lu() ;
	cout << "Matching matrix :" << endl ;
	cout << matching << endl ;
	
	// Second member for the matching matrix :
	Tbl sec_member (4) ;
	sec_member.annule_hard() ;
	// Boundary on the left :
	for (int i=0 ; i<nr ; i++)
	     sec_member.set(0) += (i%2==0) ? -sp_left(i) : sp_left(i) ;
	// Matching :
	for (int i=0 ; i<nr ; i++)
	     sec_member.set(1) += (i%2==0) ? -sp_left(i)-sp_right(i)  : -sp_left(i)+sp_right(i) ;
	// Matching of the derivative :
	for (int i=0 ; i<nr ; i++)
	     sec_member.set(2) += (i%2==0) ? i*i*(sp_left(i)+sp_right(i))  : i*i*(sp_left(i)-sp_right(i)) ;
	//Boundary on the right :
	for (int i=0 ; i<nr ; i++)
	     sec_member.set(3) -= sp_right(i) ;
	     
	cout << "Second member of the matching system: " << endl ;
	cout << sec_member << endl ;
	
	Tbl coef_sh (matching.inverse(sec_member)) ;
	cout << "Coefficients of the homogeneous solutions :" << endl ;
	cout << coef_sh << endl ;
	 
	// Coefficients of the result :
	Tbl res_left (sp_left + coef_sh(0) * coef_sh_moins_left + coef_sh(1) * coef_sh_plus_left) ;
	Tbl res_right (sp_right + coef_sh(2) * coef_sh_plus_right + coef_sh(3) * coef_sh_moins_right) ;
	
	cout << "Coefficients of the solution : " << endl ;
	cout << res_left << endl ;
	cout << res_right << endl ;
	
	// Output in a file for plotting purposes :
	char name_out[30] ;
	sprintf(name_out, "plot_multi_homogeneous_%i.dat", nr) ;
	ofstream fiche (name_out) ;
	
	int resolution = 200 ;
	double x=-1 ;
	double step = 2./resolution ;
	
	// We will also compute the maximum difference between the solution and its numerical value :
	double error_max = 0 ;
	
	double xi ;
	for (int i=0 ;  i<resolution+1 ; i++)  {
	
	    // in which domain ?
	    if (x<0)
		xi = 2*x+1 ;
	    else
		xi = 2*x-1 ;
	    
	    // One computes the solution at the current point
	    double val_func = 0 ;
	    for (int j=0 ; j<nr ; j++)
	         if (x<0)
	              val_func += res_left(j)*cheby(j, xi) ;
		 else
		      val_func += res_right(j)*cheby(j, xi) ;
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
