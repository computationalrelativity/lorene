/*
 * Method of the class Valeur to compute an equipotential surface
 *  (outward search)
 *
 * (see file valeur.h for the documentation).
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char valeur_equipot_out_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.2  2001/01/08  10:48:27  eric
 * Version moins severe avec les discontinuites entre les domaines
 *  (on a remplace un abort() par un warning).
 *
 * Revision 2.1  2000/11/10  15:20:55  eric
 * Les 1.e-14 sont remplaces par precis0 = precis.
 *
 * Revision 2.0  2000/11/10  13:32:19  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers C++
#include <iostream.h>

// Headers C
#include <stdlib.h>
#include <math.h>

// Headers Lorene
#include "valeur.h"
#include "itbl.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"

// Local prototypes
double valeur_equipot_fonc(double, const Param&) ;

//****************************************************************************
 
void Valeur::equipot_outward(double uu0, int nz_search, double precis, 
			     int nitermax, int& niter, Itbl& l_iso, 
			     Tbl& xi_iso) const {

    
    int nz = mg->get_nzone() ; 
    int nt = mg->get_nt(0) ; 
    int np = mg->get_np(0) ; 
    
    double precis0 = precis ; 
    
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
//    immediatement apres r_iso(phi,theta)

	int i2 ;
	l2 = nz ;
	for (int l=0; l<nz_search; l++) {	
	    int nr = mg->get_nr(l) ;
	
	    for (int i=0; i<nr; i++) {
		double uux = (*this)(l, k, j, i) ;
		if ( ( (uux < uu0) || ( fabs(uux-uu0) < precis0 ) ) &&
			(uux != __infinity) ) {
		    l2 = l ;		// on a trouve le point 
		    i2 = i ;		//
		    break ;
		}
	    }

	    if (l2 == l) break ;
	}   // fin de la boucle sur les zones

	if (l2==nz) {
	    cout << 
		"Valeur::equipot_outward: the point uu < uu0 has not been found"
		 << endl ;
	    cout << " for the phi index " << k 
	         << " and the theta index " << j <<  endl ;
	    cout << " uu0 = " << uu0 << endl ;
	    abort () ;
	}	    

// 2. Point precedent (l2,i3) 
	int i3 ;
	if (i2 != 0) { // on reste dans la meme zone
	    i3 = i2 - 1 ;
	}
	else {  
	     if (l2 != 0) { // on change de zone

		l2-- ; 
		i2 = mg->get_nr(l2) - 1 ; 

		double uux = (*this)(l2, k, j, i2) ;

		if ( ( fabs(uux-uu0) > precis0 ) ) {
		    // the field may slightly discontinuous:
	    	    cout << 
	"Valeur::equipot_outward: WARNING: potentially discontinuous field !" 
		    << endl ;
		    cout << " k,  j : " << k << " " << j << endl ;
		    cout << " uux,  uu0 : " << uux << "   " << uu0 << endl ;
		} 

		l_iso.set(k, j) = l2 ;
		xi_iso.set(k, j) = (mg->get_grille3d(l2))->x[i2] ;
		continue ;	// theta suivant

	     }
	     else {	// l2 = 0, i2 = 0
		cout << 
  "Valeur::equipot_outward: the field has some negative value at the center !" 
		<< endl ;
		cout << " k,  j : " << k << " " << j << endl ;
		cout << " uu0 : " << uu0 << endl ;
		abort () ;
	     }
	}

	l_iso.set(k, j) = l2 ;

	double x2 = (mg->get_grille3d(l2))->x[i2] ;
	double x3 = (mg->get_grille3d(l2))->x[i3] ;

	double y2 = (*this)(l2, k, j, i2) - uu0 ;

// 3. Recherche de x0 
	
	int niter0 = 0 ;
	if (fabs(y2) < precis0) {
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
    
