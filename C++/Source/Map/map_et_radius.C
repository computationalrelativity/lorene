/*
 *  Methods of the class Map_et relative to the function
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


char map_et_radius_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2013/06/05 15:10:42  j_novak
 * Suppression of FINJAC sampling in r. This Jacobi(0,2) base is now
 * available by setting colloc_r to BASE_JAC02 in the Mg3d constructor.
 *
 * Revision 1.4  2008/08/27 08:49:16  jl_cornou
 * Added R_JACO02 case
 *
 * Revision 1.3  2004/01/26 16:58:35  j_novak
 * Added initialization to avoid compiler warning.
 *
 * Revision 1.2  2003/12/19 16:21:43  j_novak
 * Shadow hunt
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.5  2000/10/20  13:14:47  eric
 * Map_et::val_lx : la valeur par defaut de precis est desormais 1e-14
 *  (et non plus 1e-15).
 *
 * Revision 1.4  2000/01/04  13:05:47  eric
 * \Corrections dans val_lx et val_lx_jk : initialisation de ftp/gtp
 *  par double(0) dans les cas ou il ne sont pas calcules.
 *
 * Revision 1.3  1999/12/17  11:03:36  eric
 * Fonctions val_r_jk et val_lx_jk operationnelles.
 *
 * Revision 1.2  1999/12/17  09:29:22  eric
 * val_lx : initialisation de niter a 0.
 *
 * Revision 1.1  1999/12/16  14:19:39  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */


// Headers Lorene
#include "map.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"

// Local prototypes
double fonc_invr_map_et_noyau(double, const Param&) ; 
double fonc_invr_map_et_coq(double, const Param&) ; 
double fonc_invr_map_et_zec(double, const Param&) ; 

			//------------------------------// 
			//	    val_r		//
			//------------------------------// 

 
double Map_et::val_r(int l, double xi, double theta, double pphi) const {

    assert( l>=0 ) ; 
    assert( l<mg->get_nzone() ) ; 
    
    double resu ; 
    double ftp = ff.val_point(l, 0, theta, pphi) ; // value of F_l(theta,phi)

    switch( mg->get_type_r(l) ) {

	case RARE: {
	    double gtp = gg.val_point(l, 0, theta, pphi) ; 
	    double xi_2 = xi * xi ; 
	    double xi_3 = xi * xi_2 ;
	    double a = xi_2 * xi_2 * (3. - 2.*xi_2) ;
	    double b =  ( 2.5  - 1.5 * xi_2 ) * xi_3 ;
	    resu = alpha[l] * ( xi + a * ftp + b * gtp ) + beta[l] ;
	    break ;
	}

	case FIN: {
	    double gtp = gg.val_point(l, 0, theta, pphi) ; 
	    double xm1 = xi - 1. ; 
	    double xp1 = xi + 1. ; 
	    double a = 0.25* xm1 * xm1 * (xi + 2.) ;
	    double b = 0.25* xp1 * xp1 * (2. - xi) ;
	    resu = alpha[l] * ( xi + a * ftp + b * gtp ) + beta[l] ;
	    break ;
	}
	
	case UNSURR: {
	    double xm1 = xi - 1. ; 
	    double a = 0.25* xm1 * xm1 * (xi + 2.) ;
	    resu = double(1) / ( alpha[l] * ( xi + a * ftp ) + beta[l] ) ;
	    break ;
	}

	default: {
	    cout << "Map_et::val_r: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }
             
    return resu ;    
}
			
			//------------------------------// 
			//	    val_lx		//
			//------------------------------// 

//-------------------------------
// Version with default precision 
//-------------------------------

void Map_et::val_lx(double rr, double theta, double pphi,  
			    int& lz, double& xi) const {

    int nitermax = 100 ;    // Maximum number of iteration in the secant method
    int niter ; 
    double precis = 1e-14 ; // Absolute precision in the secant method
    
    Param par ; 
    par.add_int(nitermax) ;	
    par.add_int_mod(niter) ; 
    par.add_double(precis) ; 

    // Call of the version with precision parameters
    
    val_lx(rr, theta, pphi, par, lz, xi) ;

}


//---------------------------------
// Version with specified precision 
//---------------------------------

void Map_et::val_lx(double rr, double theta, double pphi,
		    const Param& par, int& lz, double& xi) const {
			   
    int nz = mg->get_nzone() ;

    // Precision in the secant method :
    // ------------------------------
    
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    niter = 0 ;	    // initialisation 
    double precis = par.get_double() ; 
    
    // Particular case of r = infinity
    //--------------------------------
    if (rr == __infinity) { 
	if ( (mg->get_type_r(nz-1) == UNSURR) && 
	     (alpha[nz-1] == - beta[nz-1]) ) {
	    lz = nz - 1 ; 
	    xi = 1 ;		 
	}
	else {
	    cout.precision(16);
	    cout.setf(ios::showpoint);
	    cout << "Map_et::val_lx: the domain containing r = " << rr <<
		 " has not been found ! " 
	      << endl ;
	    abort () ;
	}
	return ; 
    }
    

    // In which domain is located r ? 
    // ----------------------------
    lz = - 1 ;
    double rmax = 0. ; 
    double ftp = double(0) ; 
    double gtp = double(0) ; 
    for (int l=0; l<nz; l++) {
		
	switch( mg->get_type_r(l) ) {

	case RARE: {
	    ftp = ff.val_point(l, 0, theta, pphi) ; // value of F_l(theta,phi)
	    gtp = gg.val_point(l, 0, theta, pphi) ; // value of G_l(theta,phi)
	    rmax = alpha[l] * ( double(1) + ftp + gtp ) + beta[l] ;
	    break ;
	}
	
	case FIN: {
	    ftp = double(0) ; 
	    gtp = gg.val_point(l, 0, theta, pphi) ; 
	    rmax = alpha[l] * ( double(1) + gtp ) + beta[l] ;
	    break ;
	}
	
	case UNSURR: {
	    ftp = double(0) ; 
	    gtp = double(0) ; 
	    rmax = double(1) / ( alpha[l] + beta[l] ) ; 
	    break ;
	}

	default: {
	    cout << "Map_et::val_lx: unknown type_r ! " << endl ;
	    abort() ;
	}	   
	}   // end of switch on type_r


	if ( rr <= rmax + 1.e-14 ) { 
	    lz = l ;
	    if (ftp == double(0)) ftp = ff.val_point(l, 0, theta, pphi) ;
	    if (gtp == double(0)) gtp = gg.val_point(l, 0, theta, pphi) ;
	    break ; 
	}	
    }		// End of loop onto the domains
    
    if (lz == -1) {		    // The domain has not been found
	cout.precision(16);
	cout.setf(ios::showpoint);
	cout << "Map_et::val_lx: the domain containing r = " << rr <<
		 " has not been found ! " 
	      << endl ;
	for (int l=0; l<nz; l++) {
	    ftp = ff.val_point(l, 0, theta, pphi) ;
	    gtp = gg.val_point(l, 0, theta, pphi) ; 
	    switch( mg->get_type_r(l) ) {
		case RARE: {
		    rmax = alpha[l] * ( double(1) + ftp + gtp ) + beta[l] ;
		    break ;
		}
		case FIN: {
		    rmax = alpha[l] * ( double(1) + gtp ) + beta[l] ;
		    break ;
		}
		case UNSURR: {
		    rmax = double(1) / ( alpha[l] + beta[l] ) ; 
		    break ;
		}
		default: {
		    cout << "Map_et::val_lx: unknown type_r ! " << endl ;
		    abort () ;
		}	   
	    }   // end of switch on type_r

	    cout << "domain: " << l << " theta = " << theta << " phi = " 
		 << pphi << " :  rmax = " << rmax << endl ; 
	}
	abort() ;
    }

    // Computation of xi
    // ----------------- 

    if ( (rr >= rmax) && ( rr <= rmax + 1.e-14) ) {
	xi = double(1) ; 
    }
    else {

    // Search of xi by the secant method
    // ---------------------------------

    Param parzerosec ;
    parzerosec.add_double(rr, 0) ; 
    parzerosec.add_double(ftp, 1) ; 
    parzerosec.add_double(gtp, 2) ; 
    parzerosec.add_double(alpha[lz], 3) ; 
    parzerosec.add_double(beta[lz], 4) ; 


    switch( mg->get_type_r(lz) ) {

	case RARE: {
	    if ( (ff.get_etat()==ETATZERO) && (gg.get_etat()==ETATZERO) ) {
		xi = ( rr - beta[lz] ) / alpha[lz]  ;
	    }
	    else {
		double xmin = 0 ; 
		double xmax = 1 ; 
		xi = zerosec(fonc_invr_map_et_noyau, parzerosec, xmin, xmax, 
			     precis, nitermax, niter) ;
	    }
	    break ;
	}
	
	case FIN: {
	    if ( (ff.get_etat()==ETATZERO) && (gg.get_etat()==ETATZERO) ) {
		xi = ( rr - beta[lz] ) / alpha[lz]  ;
	    }
	    else {
		double xmin = -1 ; 
		double xmax = 1 ; 
		xi = zerosec(fonc_invr_map_et_coq, parzerosec, xmin, xmax, 
			     precis, nitermax, niter) ;
	    }
	    break ;
	}
	
	case UNSURR: {
	    if ( (ff.get_etat()==ETATZERO) ) {
		xi = ( double(1)/rr - beta[lz] ) / alpha[lz]  ;
	    }
	    else {
		assert(ff.get_etat()==ETATQCQ) ; 
		if ( ff.c->t[lz]->get_etat() == ETATZERO) {
		    xi = ( double(1) / rr - beta[lz] ) / alpha[lz]  ;
		}
		else {
		    double xmin = -1 ; 
		    double xmax = 1 ; 
		    xi = zerosec(fonc_invr_map_et_zec, parzerosec, xmin, xmax, 
			     precis, nitermax, niter) ;
		} 
	    }
	    break ;
	}

	default: {
	    cout << "Map_et::val_lx: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }
    
    
    } // End of search by the secant method

} 




			//------------------------------// 
			//	    val_r_jk		//
			//------------------------------// 

 
double Map_et::val_r_jk(int l, double xi, int j, int k) const {

    assert( l>=0 ) ; 
    assert( l<mg->get_nzone() ) ; 
    
    double resu ; 
    double ftp = ff(l, k, j, 0) ; // value of F_l(theta_j, phi_k)

    switch( mg->get_type_r(l) ) {

	case RARE: {
	    double gtp = gg(l, k, j, 0) ; // value of G_l(theta_j, phi_k)
	    double xi_2 = xi * xi ; 
	    double xi_3 = xi * xi_2 ;
	    double a = xi_2 * xi_2 * (3. - 2.*xi_2) ;
	    double b =  ( 2.5  - 1.5 * xi_2 ) * xi_3 ;
	    resu = alpha[l] * ( xi + a * ftp + b * gtp ) + beta[l] ;
	    break ;
	}
	
	case FIN: {
	    double gtp = gg(l, k, j, 0) ; // value of G_l(theta_j, phi_k)
	    double xm1 = xi - 1. ; 
	    double xp1 = xi + 1. ; 
	    double a = 0.25* xm1 * xm1 * (xi + 2.) ;
	    double b = 0.25* xp1 * xp1 * (2. - xi) ;
	    resu = alpha[l] * ( xi + a * ftp + b * gtp ) + beta[l] ;
	    break ;
	}

	case UNSURR: {
	    double xm1 = xi - 1. ; 
	    double a = 0.25* xm1 * xm1 * (xi + 2.) ;
	    resu = double(1) / ( alpha[l] * ( xi + a * ftp ) + beta[l] ) ;
	    break ;
	}

	default: {
	    cout << "Map_et::val_r_jk: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }
             
    return resu ;    

}
			
			//------------------------------// 
			//	    val_lx_jk		//
			//------------------------------// 

void Map_et::val_lx_jk(double rr, int j, int k, const Param& par, 
			    int& lz, double& xi) const {
			   
    int nz = mg->get_nzone() ;

    // Precision in the secant method :
    // ------------------------------
    
    int nitermax = par.get_int() ; 
    int& niter = par.get_int_mod() ; 
    niter = 0 ;	    // initialisation 
    double precis = par.get_double() ; 
    
    // Particular case of r = infinity
    //--------------------------------
    if (rr == __infinity) { 
	if ( (mg->get_type_r(nz-1) == UNSURR) && 
	     (alpha[nz-1] == - beta[nz-1]) ) {
	    lz = nz - 1 ; 
	    xi = 1 ;		 
	}
	else {
	    cout.precision(16);
	    cout.setf(ios::showpoint);
	    cout << "Map_et::val_lx_jk: the domain containing r = " << rr <<
		 " has not been found ! " 
	      << endl ;
	    abort () ;
	}
	return ; 
    }
    

    // In which domain is located r ? 
    // ----------------------------
    lz = - 1 ;
    double rmax = 0 ; 
    double ftp = double(0) ; 
    double gtp = double(0) ; 
    for (int l=0; l<nz; l++) {
		
	switch( mg->get_type_r(l) ) {

	case RARE: {
	    ftp = ff(l, k, j, 0) ; // value of F_l(theta_j, phi_k)
	    gtp = gg(l, k, j, 0) ; // value of G_l(theta_j, phi_k)
	    rmax = alpha[l] * ( double(1) + ftp + gtp ) + beta[l] ;
	    break ;
	}
	
	case FIN: {
	    ftp = double(0) ; 
	    gtp = gg(l, k, j, 0) ; // value of G_l(theta_j, phi_k)
	    rmax = alpha[l] * ( double(1) + gtp ) + beta[l] ;
	    break ;
	}
	
	case UNSURR: {
	    ftp = double(0) ; 
	    gtp = double(0) ; 
	    rmax = double(1) / ( alpha[l] + beta[l] ) ; 
	    break ;
	}

	default: {
	    cout << "Map_et::val_lx_jk: unknown type_r ! " << endl ;
	    abort() ;
	}	   
	}   // end of switch on type_r


	if ( rr <= rmax + 1.e-14 ) { 
	    lz = l ;
	    if (ftp == double(0)) ftp = ff(l, k, j, 0) ;
	    if (gtp == double(0)) gtp = gg(l, k, j, 0) ;
	    break ; 
	}	
    }		// End of loop onto the domains
    
    if (lz == -1) {		    // The domain has not been found
	cout.precision(16);
	cout.setf(ios::showpoint);
	cout << "Map_et::val_lx_jk: the domain containing r = " << rr <<
		 " has not been found ! " 
	      << endl ;
	for (int l=0; l<nz; l++) {
	    ftp = ff(l, k, j, 0) ;
	    gtp = gg(l, k, j, 0) ; 
	    switch( mg->get_type_r(l) ) {
		case RARE: {
		    rmax = alpha[l] * ( double(1) + ftp + gtp ) + beta[l] ;
		    break ;
		}
		case FIN: {
		    rmax = alpha[l] * ( double(1) + gtp ) + beta[l] ;
		    break ;
		}
		case UNSURR: {
		    rmax = double(1) / ( alpha[l] + beta[l] ) ; 
		    break ;
		}
		default: {
		    cout << "Map_et::val_lx_jk: unknown type_r ! " << endl ;
		    abort () ;
		}	   
	    }   // end of switch on type_r

	    cout << "domain: " << l << "   j = " << j << "   k = " << k
		 << " :  rmax = " << rmax << endl ; 
	}
	abort() ;
    }

    // Computation of xi
    // ----------------- 

    if ( (rr >= rmax) && ( rr <= rmax + 1.e-14) ) {
	xi = double(1) ; 
    }
    else {

    // Search of xi by the secant method
    // ---------------------------------

    Param parzerosec ;
    parzerosec.add_double(rr, 0) ; 
    parzerosec.add_double(ftp, 1) ; 
    parzerosec.add_double(gtp, 2) ; 
    parzerosec.add_double(alpha[lz], 3) ; 
    parzerosec.add_double(beta[lz], 4) ; 


    switch( mg->get_type_r(lz) ) {

	case RARE: {
	    if ( (ff.get_etat()==ETATZERO) && (gg.get_etat()==ETATZERO) ) {
		xi = ( rr - beta[lz] ) / alpha[lz]  ;
	    }
	    else {
		double xmin = 0 ; 
		double xmax = 1 ; 
		xi = zerosec(fonc_invr_map_et_noyau, parzerosec, xmin, xmax, 
			     precis, nitermax, niter) ;
	    }
	    break ;
	}
	
	case FIN: {
	    if ( (ff.get_etat()==ETATZERO) && (gg.get_etat()==ETATZERO) ) {
		xi = ( rr - beta[lz] ) / alpha[lz]  ;
	    }
	    else {
		double xmin = -1 ; 
		double xmax = 1 ; 
		xi = zerosec(fonc_invr_map_et_coq, parzerosec, xmin, xmax, 
			     precis, nitermax, niter) ;
	    }
	    break ;
	}
	
	case UNSURR: {
	    if ( (ff.get_etat()==ETATZERO) ) {
		xi = ( double(1)/rr - beta[lz] ) / alpha[lz]  ;
	    }
	    else {
		assert(ff.get_etat()==ETATQCQ) ; 
		if ( ff.c->t[lz]->get_etat() == ETATZERO) {
		    xi = ( double(1) / rr - beta[lz] ) / alpha[lz]  ;
		}
		else {
		    double xmin = -1 ; 
		    double xmax = 1 ; 
		    xi = zerosec(fonc_invr_map_et_zec, parzerosec, xmin, xmax, 
			     precis, nitermax, niter) ;
		} 
	    }
	    break ;
	}

	default: {
	    cout << "Map_et::val_lx_jk: unknown type_r ! " << endl ;
	    abort () ;
	}	   
    }
    
    
    } // End of search by the secant method
} 


//=============================================================================//
//				Auxilary functions			       //
//=============================================================================//
	   
double fonc_invr_map_et_noyau(double x, const Param& par) {

    double r = par.get_double(0) ;
    double f = par.get_double(1) ;
    double g = par.get_double(2) ;
    double alp = par.get_double(3) ;
    double bet = par.get_double(4) ;
    double x_2 = x * x ; 
    double x_3 = x_2 * x ;
    double a = x_2 * x_2 * (3. - 2.*x_2) ;
    double b =  ( 2.5  - 1.5 * x_2 ) * x_3 ;
    
    return alp * ( x + a * f + b * g ) + bet - r ;
 
}

//****************************************************************************

double fonc_invr_map_et_coq(double x, const Param& par) {

    double r = par.get_double(0) ;
    double f = par.get_double(1) ;
    double g = par.get_double(2) ;
    double alp = par.get_double(3) ;
    double bet = par.get_double(4) ;
    double xm1 = x - 1. ;
    double xp1 = x + 1. ;
    double a = 0.25* xm1 * xm1 * (x + 2.) ;
    double b = 0.25* xp1 * xp1 * (2. - x) ;
    
    return alp * ( x + a * f + b * g ) + bet - r ;
 
}
	   
//****************************************************************************

double fonc_invr_map_et_zec(double x, const Param& par) {

    double r = par.get_double(0) ;
    double f = par.get_double(1) ;
    double alp = par.get_double(3) ;
    double bet = par.get_double(4) ;
    double xm1 = x - 1. ;
    double a = 0.25* xm1 * xm1 * (x + 2.) ;
    
    return alp * ( x + a * f ) + bet - double(1) / r ;
 
}

