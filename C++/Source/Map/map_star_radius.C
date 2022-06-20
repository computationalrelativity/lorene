/*
 *  Methods of the class Map_star relative to the function
 *	    r = R_l(xi, theta', phi')
 */

/*
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








#include <cmath>

// Headers Lorene
#include "map.h"

			//------------------------------// 
			//	          val_r	         	//
			//------------------------------// 

 
namespace Lorene {
double Map_star::val_r(int l, double xi, double theta, double pphi) const {

    assert( l>=0 ) ; 
    assert( l<mg->get_nzone() ) ; 
    
    double resu ; 

    switch( mg->get_type_r(l) ) {

	case RARE: {
	    resu = alpha.val_point(l, xi, theta, pphi) * xi ;
	    break ;
	}
    
    case FIN: {
	    resu = alpha.val_point(l, xi, theta, pphi)*xi + beta.val_point(l, xi, theta, pphi) ;
	    break ;
	}


	case UNSURR: {
	    cout << "Map_star::val_r : Compactified domain not allowed !" << endl;
        abort() ;
	    break ;
	}

	default: {
	    cout << "Map_star::val_r: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }
             
    return resu ;    
}
			
			//------------------------------// 
			//	          val_lx	     	//
			//------------------------------// 

void Map_star::val_lx(double rr, double theta, double pphi, int& lz, double& xi) const {
			   
    // In which domain is located r ? 
    // ----------------------------
    int nz = mg->get_nzone() ;
    lz = - 1 ;
    
    for (int l=0; l<nz; l++) {
	double rmin,rmax ; 

	if (mg->get_type_r(l) == RARE){
        rmin = 0. ;
        rmax = alpha.val_point(l,1.,theta,pphi) ;
    }
	if (mg->get_type_r(l) == FIN) {
		// cout << "BONJOUR" << endl;
	    rmin = beta.val_point(l,1.,theta,pphi) - alpha.val_point(l,1,theta,pphi) ;
		rmax = beta.val_point(l,1.,theta,pphi) + alpha.val_point(l,1,theta,pphi) ;
	}
	// if (mg->get_type_r(l) == UNSURR) {
	//     rmin = double(1)/rmin ; 
	//     rmax = double(1)/rmax ; 
	// }		
	if ((rr - rmin >= -1.e-14*fabs(rmin)) && ( rr <= rmax )) { 
	    lz = l ;
	    break ; 
	}	
    }		// fin de la boucle sur les zones
    
    if (lz == -1) {		    // On n'a pas trouve la zone 
	cout.precision(16);
	cout.setf(ios::showpoint);
	cout << "Map_star::val_lx: the domain containing r = " << rr <<
		 " has not been found ! " 
	      << endl ;
	// for (int l=0; l<nz; l++) {
	//     double rmin = -alpha[l] + beta[l] ;	
	//     if (mg->get_type_r(l) == UNSURR) rmin = double(1)/rmin ; 
	//     if (mg->get_type_r(l) == RARE) rmin = 0. ;
	//     cout << "domain " << l << " :  r_min = " << rmin ; 
	//     double rmax = alpha[l] + beta[l] ;	
	//     if (mg->get_type_r(l) == UNSURR) rmax = double(1)/rmax ; 
	//     cout << " :  r_max = " << rmax << endl ; 
	// }
	abort () ;
    }

    // Computation of xi
    // ----------------- 

    switch( mg->get_type_r(lz) ) {

	case RARE: {
	    xi = rr / alpha.val_point(lz,1.,theta,pphi)  ;
	    break ;
	}
	
	case FIN: {
	    xi = (rr-beta.val_point(lz,1.,theta,pphi))/alpha.val_point(lz,1.,theta,pphi) ;
	    break ;
	}


	case UNSURR: {
	    cout << "Map_star::val_lx : Compactified domain not allowed !" << endl;
        abort() ;
	    break ;
	}

	default: {
	    cout << "Map_star::val_lx: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }

} 


void Map_star::val_lx(double rr, double theta, double phi, const Param&,  
			    int& lz, double& xi) const {

    val_lx(rr, theta, phi, lz, xi) ;

}


			//------------------------------// 
			//	    val_r_jk		        //
			//------------------------------// 

 
double Map_star::val_r_jk(int l, double xi, int j, int k) const {

    return val_r(l, xi, (+tet)(l,k,j,0), (+phi)(l,k,j,0)) ; 
    
}
			
			//------------------------------// 
			//	    val_lx_jk		        //
			//------------------------------// 

void Map_star::val_lx_jk(double rr, int j, int k, const Param& par,
			    int& l, double& xi) const {
			   
    val_lx(rr, (+tet)(l,k,j,0), (+phi)(l,k,j,0), l, xi) ; 

} 


}