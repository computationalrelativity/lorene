/*
 * Computation of Legendre Gauss-Lobatto nodes and weights
 */
 
/*
 *   Copyright (c) 2005 Philippe Grandclément
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/11/14 01:56:59  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

#include <iostream>

using namespace std ;

#include "math.h"

void legendre_poly_der(int n, double& poly, double& pder, double& polym1, 
        double& pderm1, double& polym2, double& pderm2, double x) {
		
	
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


void legendre_nodes_weight_GL(int n, double* coloc, double* weight, double prec,
                              int itemax) {
     
     double x_plus = 1 ;
     double x_moins = -1 ;
     double p_plus, dp_plus, p_plus_m1, dp_plus_m1, p_plus_m2, dp_plus_m2 ;
     double p_moins, dp_moins, p_moins_m1, dp_moins_m1, p_moins_m2, dp_moins_m2 ;
     double p, dp, p_m1, dp_m1, p_m2, dp_m2 ;
     
     legendre_poly_der(n, p_plus, dp_plus, p_plus_m1, dp_plus_m1, p_plus_m2, 
                       dp_plus_m2, x_plus) ;
     legendre_poly_der(n, p_moins, dp_moins, p_moins_m1, dp_moins_m1, 
                       p_moins_m2, dp_moins_m2, x_moins) ;
     
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
	       legendre_poly_der(n, p, dp, p_m1, dp_m1, p_m2, dp_m2, x) ;
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
	    cout << "legendre_poly_der: too many iterations !..." << endl ;
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
          legendre_poly_der(n-1, p, dp, p_m1, dp_m1, p_m2, dp_m2, coloc[i]) ;
	  weight[i] = 2./(n-1)/n/p/p ;
    }
}
	

