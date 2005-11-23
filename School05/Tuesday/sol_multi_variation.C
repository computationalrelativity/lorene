#include "headcpp.h"
#include "math.h"

#include "tbl.h"
#include "matrice.h"
#include "leg.h"

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
	for (int n=0 ; n<size ; n++) {
	     for (int p=n+1 ; p<size ; p++)
		if ((p+n)%2 == 1)
			res.set(n) += so(p) ;
		res.set(n) *= (2*n+1) ;
	}
	return res ;
}


int main () {

	int nr ;
	cout << "Please enter the number of coefficients :" << endl ;
	cin >> nr ;
	
	// Collocation and weights
	Tbl coloc (nr) ;
	coloc.set_etat_qcq() ;
	Tbl poids (nr) ;
	poids.set_etat_qcq() ;
	coloc_poids_leg(nr, coloc, poids) ;

	// Computation of the matrix D :
	Matrice ope_D (nr, nr) ;
	ope_D.set_etat_qcq() ;
	
	Tbl val (nr) ;
	for (int i=0 ; i<nr ; i++) {
		val.annule_hard() ;
		val.set(i) = 1 ;
		Tbl coef (coef_leg(val)) ;
		Tbl der (ope_der(coef)) ;
		
		// Value of derivative at point x_j :
		for (int j=0 ; j<nr ; j++) {
			double sum = 0 ;
			for (int k=0 ; k<nr ; k++)
		    		sum += der(k)*leg(k, coloc(j)) ;
 			ope_D.set(j,i) = 2*sum ;
		} 
	}
	
	cout << "Derivative operator D : " << endl ;
	cout << ope_D << endl ;
	
	// Matrice of the system :
	Matrice systeme(2*nr-1, 2*nr-1) ;
	systeme.annule_hard() ;
	
	// Left boundary condition : 
	systeme.set(0,0) = 1 ;
	// Equation inside left domain and at the interface:
	for (int n=1 ; n<nr ; n++) {
	       for (int j=0 ; j<nr ; j++)
	            for (int i=0 ; i<nr ; i++)
		         systeme.set(n,j) += ope_D(i,j)*ope_D(i,n)*poids(i) ;
		systeme.set(n,n) += 4*poids(n) ;
	}
	// equation at the interface and inside domain 1 :
	for (int n=0 ; n<nr-1 ; n++) {
	       for (int j=0 ; j<nr ; j++)
	            for (int i=0 ; i<nr ; i++)
		         systeme.set(n+nr-1,j+nr-1) += ope_D(i,j)*ope_D(i,n)*poids(i) ;
		systeme.set(n+nr-1,n+nr-1) += 4*poids(n) ;
	}
	// Boundary on the right : 
	systeme.set(2*nr-2, 2*nr-2) = 1 ;

	cout << "Matrix of the system: " << endl ;
	cout << systeme << endl ;
	systeme.set_lu() ;
	
	// Construction of the second member :
	Tbl sec_member (2*nr-1) ;
	sec_member.annule_hard() ;
	// Left boundary condition :
	sec_member.set(0) = 0 ;
	// domain 1 ;
	for (int i=1 ; i<nr ; i++)
	    sec_member.set(i) = f_left(coloc(i))*poids(i) ;
	// domaine 2 ;
	for (int i=0 ; i<nr-1 ; i++)
	    sec_member.set(i+nr-1) += f_right(coloc(i))*poids(i) ;
	// right boundary condition :
	sec_member.set(2*nr-2) = 0 ;
	
	cout << "Second member of the system" << endl ;
	cout << sec_member << endl ;
	
	// Inversion of the system :
	Tbl inv (systeme.inverse(sec_member)) ;
	// Solution in both domain, at collocations points :
	Tbl sol_left (nr) ;
	sol_left.set_etat_qcq() ;
	for (int i=0 ; i<nr ; i++)
	     sol_left.set(i) = inv(i) ;
	Tbl sol_right (nr) ;
	sol_right.set_etat_qcq() ;
	for (int i=0 ; i<nr ; i++)
	     sol_right.set(i) = inv(nr-1+i) ;
	
	cout << "Solutions in both domains, at collocation points : "<< endl ;
	cout << sol_left << endl ;
	cout << sol_right << endl ;
	
	// Compute the coefficients :
	Tbl coef_left (coef_leg(sol_left)) ;
	Tbl coef_right (coef_leg(sol_right)) ;
	
	// Outputs :
	// Output in a file for plotting purposes :
	char name_out[30] ;
	sprintf(name_out, "plot_variational_%i.dat", nr) ;
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
	              val_func += coef_left(j)*leg(j, xi) ;
		 else
		      val_func += coef_right(j)*leg(j, xi) ;
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
