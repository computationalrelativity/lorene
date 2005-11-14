#include "headcpp.h"
#include "math.h"

#include "tbl.h"
#include "matrice.h"
#include "leg.h"

void legendre (int n, double& poly, double& pder, double& polym1, double& pderm1, 
		double& polym2, double& pderm2, double x) {
		
	
	if (n==0) {
	     poly = 1 ;
	     pder = 0 ;
	     }
	else 
	    if (n==1) {
	         polym1 = 1 ;
		 pderm1 = 0 ;
		 poly = x ;
		 pder = 1 ;
		 }
	else {
	     polym1 = 1 ;
	     pderm1 = 0 ;
	     poly = x ;
	     pder = 1 ;
	     for (int i=1 ; i<n ; i++) {
	         polym2 = polym1 ;
		 pderm2 = pderm1 ;
		 polym1 = poly ;
		 pderm1 = pder ;
		 poly = ((2*i+1)*x*polym1 - i*polym2)/(i+1) ;
		 pder = ((2*i+1)*polym1+(2*i+1)*x*pderm1-i*pderm2)/(i+1) ;
		}
	}
}

void gauss_lobato_legendre (int n, double* coloc, double* weight, double prec = 1e-12, int itemax = 100) {
     
     double x_plus = 1 ;
     double x_moins = -1 ;
     double p_plus, dp_plus, p_plus_m1, dp_plus_m1, p_plus_m2, dp_plus_m2 ;
     double p_moins, dp_moins, p_moins_m1, dp_moins_m1, p_moins_m2, dp_moins_m2 ;
     double p, dp, p_m1, dp_m1, p_m2, dp_m2 ;
     
     legendre (n, p_plus, dp_plus, p_plus_m1, dp_plus_m1, p_plus_m2, dp_plus_m2, x_plus) ;
     legendre (n, p_moins, dp_moins, p_moins_m1, dp_moins_m1, p_moins_m2, dp_moins_m2, x_moins) ;
     
     double det = p_plus_m1*p_moins_m2 - p_moins_m1*p_plus_m2 ;
     double r_plus = -p_plus ;
     double r_moins = -p_moins ;
     double a = (r_plus*p_moins_m2 - r_moins*p_plus_m2)/det ;
     double b = (r_moins*p_plus_m1 - r_plus*p_moins_m1)/det ;
     
     coloc[n-1] = 1 ;
     double dth = M_PI/(2*n+1) ;
     double cd = cos (2*dth) ;
     double sd = sin (2*dth) ;
     double cs = cos(dth) ;
     double ss = sin(dth) ;
     
     int borne_sup = (n%2==0) ? 
          n/2 : (n+1)/2  ;
     
     for (int j=1 ; j<borne_sup ; j++) {
          double x = cs ;
	  bool loop = true ;
	  int ite = 0 ;
	  while (loop) {
	       legendre (n, p, dp, p_m1, dp_m1, p_m2, dp_m2, x) ;
	       double poly = p + a*p_m1 + b*p_m2 ;
	       double pder = dp + a * dp_m1 + b*dp_m2 ;
	       double sum = 0 ;
	       for (int i=0 ; i<j ; i++)
	            sum += 1./(x-coloc[n-i-1]) ;
		   
	       double increm = -poly/(pder-sum*poly) ;
	       
	       x += increm ;
	       ite ++ ;
	       if ((fabs(increm) < prec) || (ite >itemax))
	            loop = false ;
	}
	if (ite > itemax) {
	    cout << "Too many iterations..." << endl ;
	    abort() ;
	}
	coloc[n-j-1] = x ;
	double auxi = cs*cd-ss*sd ;
	ss = cs*sd+ss*cd ;
	cs = auxi ;
    }
    if  (n%2==1)
        coloc[(n-1)/2] = 0 ; 
    // Copy of the symetric ones :
    for (int i=0 ; i<borne_sup ; i++)
         coloc[i] = - coloc[n-i-1] ;
      
    for (int i=0 ; i<n ; i++) {
          legendre (n-1, p, dp, p_m1, dp_m1, p_m2, dp_m2, coloc[i]) ;
	  weight[i] = 2./(n-1)/n/p/p ;
    }
}
	
// Legendre polynomial
double leg (int n, double x) {

    double p, dp, p1, dp1, p2, dp2 ;
    legendre (n, p, dp, p1, dp1, p2, dp2, x) ;
    return p ;
}

// Norme
double norme_leg (int i) {
	return  2./(2*i+1);
}

// Colocation point and weights
void coloc_poids_leg (int n, Tbl& coloc, Tbl& poids) {
     
    double* auxi_coloc = new double[n] ;
    double* auxi_poids = new double[n] ;
    gauss_lobato_legendre (n, auxi_coloc, auxi_poids) ;
    
    for (int i=0 ; i<n ; i++) {
         coloc.set(i) = auxi_coloc[i] ;
	 poids.set(i) = auxi_poids[i] ;
    }
    
    delete [] auxi_coloc ;
    delete [] auxi_poids ;
}

// Normalisation factor :
Tbl gamma_leg (int n) {
	Tbl coloc (n) ;
	coloc.set_etat_qcq() ;
	Tbl poids (n) ;
	poids.set_etat_qcq() ;
        coloc_poids_leg(n, coloc, poids) ;
	
	Tbl res (n) ;
	res.annule_hard() ;
	for (int i=0 ; i<n ; i++)
	    for (int j=0 ; j<n ; j++)
 	   	res.set(i) += leg(i,coloc(j)) *  leg(i,coloc(j)) * poids(j) ;
	return res ;
}

// Coefficients :
Tbl coef_leg (const Tbl& so) {
	
	// Verification of size :
	assert (so.get_ndim() == 1) ;
	int n = so.get_dim(0) ;
	
	// Auxiliary quantities :
	Tbl colocation (n) ;
	colocation.set_etat_qcq() ;
	Tbl poids (n) ;
	poids.set_etat_qcq() ;
	
	coloc_poids_leg(n, colocation, poids) ;
	Tbl norme (gamma_leg(n)) ;
	
	// The result is set to zero
	Tbl res (n) ;
	res.annule_hard() ;

	for (int i=0 ; i<n ; i++)
	    for (int j=0 ; j<n ; j++) 
	     	res.set(i) += so(j)*leg(i, colocation(j))*poids(j)/norme(i) ;
	return res ;
}
