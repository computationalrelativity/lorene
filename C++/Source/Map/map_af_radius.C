/*
 *  Methods of the class Map_af relative to the function
 *	    r = R_l(xi, theta', phi')
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char map_af_radius_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.5  1999/12/16  14:19:08  eric
 * Introduction de l'argument const Param& par dans val_lx et val_lx_jk.
 * (en remplacement de l'argument Tbl& param).
 *
 * Revision 1.4  1999/12/07  14:51:37  eric
 * val_r_kj --> val_r_jk
 * val_lx_kj -->val_lx_jk
 * Changement ordre des arguments val_r, val_lx
 *
 * Revision 1.3  1999/12/06  16:47:21  eric
 * Surcharge de val_lx avec la version sans param.
 *
 * Revision 1.2  1999/12/06  15:34:06  eric
 * Ajout des fonctions val_r_kj et val_lx_kj.
 *
 * Revision 1.1  1999/12/06  13:12:16  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "map.h"

			//------------------------------// 
			//	    val_r		//
			//------------------------------// 

 
double Map_af::val_r(int l, double xi, double, double) const {

    assert( l>=0 ) ; 
    assert( l<mg->get_nzone() ) ; 
    
    double resu ; 

    switch( mg->get_type_r(l) ) {

	case FIN: case RARE: {
	    resu = alpha[l] * xi + beta[l] ;
	    break ;
	}
	
	case UNSURR: {
	    resu = double(1) / ( alpha[l] * xi + beta[l] ) ;
	    break ;
	}

	default: {
	    cout << "Map_af::val_r: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }
             
    return resu ;    
}
			
			//------------------------------// 
			//	    val_lx		//
			//------------------------------// 

void Map_af::val_lx(double rr, double, double, int& lz, double& xi) const {
			   
    // In which domain is located r ? 
    // ----------------------------
    int nz = mg->get_nzone() ;
    lz = - 1 ;
    
    for (int l=0; l<nz; l++) {
        
	double rmax = alpha[l] + beta[l] ;	
	if (mg->get_type_r(l) == UNSURR) rmax = double(1)/rmax ; 
		
	if ( rr <= rmax ) { 
	    lz = l ;
	    break ; 
	}	
    }		// fin de la boucle sur les zones
    
    if (lz == -1) {		    // On n'a pas trouve la zone 
	cout.precision(16);
	cout.setf(ios::showpoint);
	cout << "Map_af::val_lx: the domain containing r = " << rr <<
		 " has not been found ! " 
	      << endl ;
	for (int l=0; l<nz; l++) {
	    double rmax = alpha[l] + beta[l] ;	
	    if (mg->get_type_r(l) == UNSURR) rmax = double(1)/rmax ; 
	    cout << "domain " << l << " :  r_max = " << rmax << endl ; 
	}
	abort () ;
    }

    // Computation of xi
    // ----------------- 

    switch( mg->get_type_r(lz) ) {

	case FIN: case RARE: {
	    xi = ( rr - beta[lz] ) / alpha[lz]  ;
	    break ;
	}
	
	case UNSURR: {
	    xi = ( double(1)/rr - beta[lz] ) / alpha[lz]  ;
	    break ;
	}

	default: {
	    cout << "Map_af::val_lx: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }

} 


void Map_af::val_lx(double rr, double, double, const Param&,  
			    int& lz, double& xi) const {

    val_lx(rr, 0., 0., lz, xi) ;

}


			//------------------------------// 
			//	    val_r_jk		//
			//------------------------------// 

 
double Map_af::val_r_jk(int l, double xi, int, int) const {

    return val_r(l, xi, 0., 0.) ; 
    
}
			
			//------------------------------// 
			//	    val_lx_jk		//
			//------------------------------// 

void Map_af::val_lx_jk(double rr, int, int, const Param& par, 
			    int& l, double& xi) const {
			   
    val_lx(rr, 0., 0., par, l, xi) ; 

} 


