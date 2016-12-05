/*
 *  Methods of the class Map_log relative to the function
 *	    r = R_l(xi, theta', phi')
 */

/*
 *   Copyright (c) 2004 Philippe Grandclement
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
 * Revision 1.5  2016/12/05 16:17:58  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:06  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:13  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2004/06/22 12:20:17  j_novak
 * *** empty log message ***
 *
 * Revision 1.1  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * 
 * $Header$
 *
 */

#include <cmath>

// Headers Lorene
#include "map.h"

			//------------------------------// 
			//	    val_r		//
			//------------------------------// 

 
namespace Lorene {
double Map_log::val_r(int l, double xi, double, double) const {

    assert( l>=0 ) ; 
    assert( l<mg->get_nzone() ) ; 
    
    double resu ; 

    switch (type_var(l)) {
    case AFFINE : {
      
      switch( mg->get_type_r(l) ) {
      case FIN: case RARE: {
	resu = alpha(l) * xi + beta(l) ;
	break ;
      }
	
      case UNSURR: {
	resu = double(1) / ( alpha(l) * xi + beta(l) ) ;
	break ;
      }
	
      default: {
	cout << "Map_log::val_r: unknown type_r ! " << endl ;
	abort () ;
      }	   
      }
      break ;
    }

    case LOG : { 
      switch( mg->get_type_r(l) ) {
      case FIN: {
	resu = exp(alpha(l) * xi + beta(l)) ;
	break ;
      }
	
      default: {
	cout << "Map_log::val_r: unknown type_r ! " << endl ;
	abort () ;
      }	   
      }
      break ;
    }
      
    default: {
      cout << "Map_log::val_r: unknown type_r ! " << endl ;
      abort () ;
    }
    }
	         
    return resu ;    
}
			
			//------------------------------// 
			//	    val_lx		//
			//------------------------------// 

void Map_log::val_lx(double rr, double, double, int& lz, double& xi) const {
			   
    // In which domain is located r ? 
    // ----------------------------
    int nz = mg->get_nzone() ;
    lz = - 1 ;
    
    for (int l=0; l<nz; l++) {
        
	double rmax = 0; 
	switch (type_var(l)) {
	case AFFINE : {
	  rmax = alpha(l) + beta(l) ;
	  break ;
	}
	case LOG : {
	  rmax = exp(alpha(l) + beta(l)) ;
	  break ;
	}
	default : {
	  cout << "Case unknown in Map_log::val_lx" << endl ;
	  break ;
	}
	}
	
	if (mg->get_type_r(l) == UNSURR) rmax = double(1)/rmax ; 
		
	if ( rr <= rmax ) { 
	    lz = l ;
	    break ; 
	}	
    }		// fin de la boucle sur les zones
    
    if (lz == -1) {		    // On n'a pas trouve la zone 
	cout.precision(16);
	cout.setf(ios::showpoint);
	cout << "Map_log::val_lx: the domain containing r = " << rr <<
	  " has not been found ! " 
	     << endl ;
	abort () ;
    }

    // Computation of xi
    // ----------------- 

    switch (type_var(lz)) {
    case AFFINE: {
      switch( mg->get_type_r(lz) ) {
      case FIN: case RARE: {
	xi = ( rr - beta(lz) ) / alpha(lz)  ;
	break ;
      }
	
      case UNSURR: {
	xi = ( double(1)/rr - beta(lz) ) / alpha(lz)  ;
	break ;
      }
	
      default: {
	cout << "Map_log::val_lx: unknown type_r ! " << endl ;
	abort () ;
      }	   
      }
      break ;
    }
    case LOG :{
      switch( mg->get_type_r(lz) ) {
      case FIN: {
	xi = ( log(rr) - beta(lz) ) / alpha(lz)  ;
	break ;
      }
      default: {
	cout << "Map_log::val_lx: unknown type_r ! " << endl ;
	abort () ;
      }	   	
      }
      break ;
    }
    default : {
      cout << "Map_log::val_lx: unknown type_r ! " << endl ;
      abort () ;
    }
    }	
}


void Map_log::val_lx(double rr, double, double, const Param&,  
			    int& lz, double& xi) const {

    val_lx(rr, 0., 0., lz, xi) ;

}


			//------------------------------// 
			//	    val_r_jk		//
			//------------------------------// 

 
double Map_log::val_r_jk(int l, double xi, int, int) const {

    return val_r(l, xi, 0., 0.) ; 
    
}
			
			//------------------------------// 
			//	    val_lx_jk		//
			//------------------------------// 

void Map_log::val_lx_jk(double rr, int, int, const Param& par, 
			    int& l, double& xi) const {
			   
    val_lx(rr, 0., 0., par, l, xi) ; 

} 


}
