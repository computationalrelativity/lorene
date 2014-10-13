/*
 * Method of the class Valeur to compute an equipotential surface.
 *
 * (see file valeur.h for the documentation).
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


char valeur_equipot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2014/10/13 08:53:50  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:23  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2003/12/19 15:05:56  j_novak
 * Added some initializations
 *
 * Revision 1.2  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/01/28  17:06:37  eric
 * Changement du cas l2<nz_search-1.
 *
 * Revision 1.2  2000/01/03  15:07:50  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/29  13:12:08  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "valeur.h"
#include "itbl.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"

// Local prototypes
namespace Lorene {
double valeur_equipot_fonc(double, const Param&) ;

//****************************************************************************
 
void Valeur::equipot(double uu0, int nz_search, double precis, int nitermax,
		     int& niter, Itbl& l_iso, Tbl& xi_iso) const {

    
    int nz = mg->get_nzone() ; 
    int nt = mg->get_nt(0) ; 
    int np = mg->get_np(0) ; 
    
    // Protections 
    // -----------
    assert(etat == ETATQCQ) ; 
    assert(nz_search > 0) ; 
    assert(nz_search <= nz) ; 
    for (int l=1; l<nz_search; l++) {
	assert( mg->get_nt(l) == nt ) ; 
	assert( mg->get_np(l) == np ) ; 
    }
    
    coef() ;	// the spectral coefficients are required
    
    l_iso.set_etat_qcq() ; 
    xi_iso.set_etat_qcq() ; 

    Param parf ;
    int j, k, l2 ; 
    parf.add_int_mod(j, 0) ; 
    parf.add_int_mod(k, 1) ; 
    parf.add_int_mod(l2, 2) ; 
    parf.add_double(uu0) ; 
    parf.add_mtbl_cf(*c_cf) ; 
    
    
    // Loop of phi and theta
    // ---------------------
    
    niter = 0 ; 
    
    for (k=0; k<np; k++) {

        for (j=0; j<nt; j++) {


// 1. Recherche du point de collocation en r (l2,i2) egal a ou situe 
//    immediatement avant r_iso(phi,theta)

	int i2 = 0;
	l2 = nz ;
	for (int l=nz_search-1; l>= 0; l--) {	// On part de la zone la plus extreme
	    int nr = mg->get_nr(l) ;
	
	    for (int i=nr-1; i >= 0; i--) {
		double uux = (*this)(l, k, j, i) ;
		if ( ( (uux > uu0) || ( fabs(uux-uu0) < 1.e-14 ) ) &&
			(uux != __infinity) ) {
		    l2 = l ;		// on a trouve le point 
		    i2 = i ;		//
		    break ;
		}
	    }

	    if (l2 == l) break ;
	}   // fin de la boucle sur les zones

	if (l2==nz) {
	    cout << "Valeur::equipot: the point uu >= uu0 has not been found" << endl ;
	    cout << " for the phi index " << k << endl ;
	    cout << " and the theta index " << j <<  endl ;
	    cout << " uu0 = " << uu0 << endl ;
	    abort () ;
	}	    

// 2. Point suivant (l2,i3) 
	int i3 = 0 ;
	if (i2 < mg->get_nr(l2) -1) { // on reste dans la meme zone
	    i3 = i2 + 1 ;
	}
	else {
	     if (l2<nz_search-1) { // on change de zone

		double uux = (*this)(l2, k, j, i2) ;
		if ( ( fabs(uux-uu0) < 1.e-14 ) ) {
		    // it is OK
	    	    cout << 
		    "Valeur::equipot: WARNING : fabs(uux-uu0) < 1.e-14" << endl ;
		    l_iso.set(k, j) = l2 ;
		    xi_iso.set(k, j) = (mg->get_grille3d(l2))->x[i2] ;
		}
		else{
	    	    cout << "Valeur::equipot: PROBLEM !!!" << endl ;
		    cout << " k,  j : " << k << " " << j << endl ;
		    cout << " uu0 : " << uu0 << endl ;
		    abort () ;
		} 
	     }
	     else {	// on est a l'extremite de la grille
		l_iso.set(k, j) = nz_search-1 ;
		xi_iso.set(k, j) = (mg->get_grille3d(l2))->x[i2] ;
		continue ;	// on passe au theta suivant
	     }
	}

	l_iso.set(k, j) = l2 ;

	double x2 = (mg->get_grille3d(l2))->x[i2] ;
	double x3 = (mg->get_grille3d(l2))->x[i3] ;

	double y2 = (*this)(l2, k, j, i2) - uu0 ;

// 3. Recherche de x0 
	
	int niter0 = 0 ;
	if (fabs(y2) < 1.e-14) {
	    xi_iso.set(k, j) = x2 ;
	}
	else {
	    xi_iso.set(k, j) = zerosec(valeur_equipot_fonc, parf, x2, x3, precis, 
				        nitermax, niter0 ) ;
	}
	
	niter = ( niter0 > niter ) ? niter0 : niter ;
	
	}	   // fin de la boucle sur theta
    }	   // fin de la boucle sur phi
        
}
    
  
		
//****************************************************************************
 
double valeur_equipot_fonc(double xi, const Param& par) {
    
    int j = par.get_int_mod(0) ; 
    int k = par.get_int_mod(1) ; 
    int l = par.get_int_mod(2) ;
    double uu0 = par.get_double() ; 
    const Mtbl_cf& cuu = par.get_mtbl_cf() ; 
    
    return cuu.val_point_jk(l, xi, j, k) - uu0 ;     
}

}
