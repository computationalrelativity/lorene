/*
 *  Arithmetics functions for the Tenseur_sym class.
 *
 *  These functions are not member functions of the Tenseur_sym class.
 *
 *  (see file tenseur.h for documentation).
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char tenseur_sym_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.3  2000/02/09  19:30:36  eric
 * MODIF IMPORTANTE: la triade de decomposition est desormais passee en
 * argument des constructeurs.
 *
 * Revision 2.2  2000/02/08  19:06:40  eric
 * Les fonctions arithmetiques ne sont plus amies.
 * Modif de diverses operations (notament division avec double)
 * Ajout de nouvelles operations (par ex. Tenseur + double, etc...)
 *
 * Revision 2.1  2000/01/11  11:15:00  eric
 * Gestion de la base vectorielle (triad).
 *
 * Revision 2.0  1999/12/02  17:18:52  phil
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
#include <assert.h>
#include <math.h>

// Headers Lorene
#include "tenseur.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

Tenseur_sym operator+(const Tenseur_sym & t) {

    return t ; 

}


Tenseur_sym operator-(const Tenseur_sym & t) {
    
   assert (t.get_etat() != ETATNONDEF) ;
    if (t.get_etat() == ETATZERO)
	return t ;
    else { 
	Tenseur_sym res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
			*(t.get_triad()) ) ; 

	res.set_etat_qcq();
	for (int i=0 ; i<t.get_n_comp() ; i++) {
	    Itbl indices (t.donne_indices(i)) ;    
	    res.set(indices) = -t(indices) ;
	    }
	return res ;
	}
}
	    

			//**********//
			// ADDITION //
			//**********//

Tenseur_sym operator+(const Tenseur_sym & t1, const Tenseur_sym & t2) {
    
    assert ((t1.get_etat() != ETATNONDEF) && (t2.get_etat() != ETATNONDEF)) ;
    assert (t1.get_valence() == t2.get_valence()) ;
    assert (t1.get_mp() == t2.get_mp()) ;
    if (t1.get_valence() != 0) {
	assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    for (int i=0 ; i<t1.get_valence() ; i++)
	assert(t1.get_type_indice(i) == t2.get_type_indice(i)) ;
   
    if (t1.get_etat() == ETATZERO)
	return t2 ;
    else if (t2.get_etat() == ETATZERO)
	    return t1 ;
	 else {
	    Tenseur_sym res(*(t1.get_mp()), t1.get_valence(), 
			    t1.get_type_indice(), *(t1.get_triad()) ) ; 

	    res.set_etat_qcq() ;
	    for (int i=0 ; i<res.get_n_comp() ; i++) {
		Itbl indices (res.donne_indices(i)) ;
		res.set(indices) = t1(indices) + t2(indices) ;
		}
	return res ;
	}
}



			//**************//
			// SOUSTRACTION //
			//**************//


Tenseur_sym operator-(const Tenseur_sym & t1, const Tenseur_sym & t2) {

    return (t1 + (-t2)) ;

}


			//****************//
			// MULTIPLICATION //
			//****************//

Tenseur_sym operator*(double x, const Tenseur_sym& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    if ( (t.get_etat() == ETATZERO) || (x == double(1)) )
	return t ;
    else {
	Tenseur_sym res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
			*(t.get_triad()) ) ; 

	if ( x == double(0) )
	    res.set_etat_zero() ;
	else {
	    res.set_etat_qcq() ;
	    for (int i=0 ; i<res.get_n_comp() ; i++) {
		Itbl indices (t.donne_indices(i)) ;
		res.set(indices) = x*t(indices) ;
		}
	    }
	    return res ; 
	}
}


Tenseur_sym operator* (const Tenseur_sym& t, double x) {
    return x * t ;
}

Tenseur_sym operator*(int m, const Tenseur_sym& t) {
    return double(m) * t ; 
}


Tenseur_sym operator* (const Tenseur_sym& t, int m) {
    return double(m) * t ;
}



			//**********//
			// DIVISION //
			//**********//

Tenseur_sym operator/ (const Tenseur_sym& t1, const Tenseur& t2) {
    
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t2.get_valence() == 0) ; // t2 doit etre un scalaire !
    assert(t1.get_mp() == t2.get_mp()) ;

    
    // Cas particuliers
    if (t2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Tenseur_sym / Tenseur !" << endl ;
	abort() ; 
    }
    if (t1.get_etat() == ETATZERO) {
    	return t1 ;
    }

    // Cas general
    
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    Tenseur_sym res(*(t1.get_mp()), t1.get_valence(), t1.get_type_indice(), 
		    *(t1.get_triad()) ) ; 

    res.set_etat_qcq() ;
    for (int i=0 ; i<res.get_n_comp() ; i++) {
	Itbl indices (res.donne_indices(i)) ;
	res.set(indices) = t1(indices) / t2() ;	    // Cmp / Cmp
    }
    return res ;

}


Tenseur_sym operator/ (const Tenseur_sym& t, double x) {

    assert (t.get_etat() != ETATNONDEF) ;
 
    if ( x == double(0) ) {
	cout << "Division by 0 in Tenseur_sym / double !" << endl ;
	abort() ;
    }

    if ( (t.get_etat() == ETATZERO) || (x == double(1)) )
	return t ;
    else {
	Tenseur_sym res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
			*(t.get_triad()) ) ; 

	res.set_etat_qcq() ;
	for (int i=0 ; i<res.get_n_comp() ; i++) {
	    Itbl indices (t.donne_indices(i)) ;
	    res.set(indices) = t(indices) / x ;	    // Cmp / double
	}
	return res ; 
    }
}



Tenseur_sym operator/ (const Tenseur_sym& t, int m) {

    return t / double(m) ; 
}


