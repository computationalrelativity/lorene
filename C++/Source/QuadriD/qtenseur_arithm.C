/*
 *  Arithmetics functions for the Qtenseur class.
 *
 *  These functions are not member functions of the Qtenseur class.
 *
 *  (see file qtenseur.h for documentation).
 *
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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

char qtenseur_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:37:13  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1  2002/09/19 09:52:43  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// Headers Lorene
#include "qtenseur.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

Qtenseur operator+(const Qtenseur & t) {

    return t ; 

}

Qtenseur operator-(const Qtenseur & t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    if (t.get_etat() == ETATZERO)
	return t ;
    else { 
	Qtenseur res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
		    t.get_triad()) ;


	res.set_etat_qcq();

	for (int i=0 ; i<res.get_n_comp() ; i++) {
	    Itbl indices (res.donne_indices(i)) ;    
	    res.set(indices) = -t(indices) ;
	    }
	return res ;
	}
}

			//**********//
			// ADDITION //
			//**********//

Qtenseur operator+(const Qtenseur & t1, const Qtenseur & t2) {
    
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
	    Qtenseur res(*(t1.get_mp()), t1.get_valence(), t1.get_type_indice(), 
			t1.get_triad() ) ;

	    res.set_etat_qcq() ;
	    for (int i=0 ; i<res.get_n_comp() ; i++) {
		Itbl indices (res.donne_indices(i)) ;
		res.set(indices) = t1(indices) + t2(indices) ;
		}
	return res ;
	}
}


Qtenseur operator+(const Qtenseur & t1, double x) {
    
    assert (t1.get_etat() != ETATNONDEF) ;
    assert (t1.get_valence() == 0) ;
    
    if (x == double(0)) {
	return t1 ;
    }
    
    Qtenseur res( *(t1.get_mp()) ) ;

    res.set_etat_qcq() ;

    res.set() = t1() + x ;	// Cmp + double

    return res ;

}


Qtenseur operator+(double x, const Qtenseur & t2) {
    
    return t2 + x ; 
    
}

Qtenseur operator+(const Qtenseur & t1, int m) {
    
    return t1 + double(m) ;

}


Qtenseur operator+(int m, const Qtenseur & t2) {
    
    return t2 + double(m) ; 
    
}



			//**************//
			// SOUSTRACTION //
			//**************//

Qtenseur operator-(const Qtenseur & t1, const Qtenseur & t2) {

    return (t1 + (-t2)) ;

}


Qtenseur operator-(const Qtenseur & t1, double x) {

    assert (t1.get_etat() != ETATNONDEF) ;
    assert (t1.get_valence() == 0) ;
    
    if (x == double(0)) {
	return t1 ;
    }
    
    Qtenseur res( *(t1.get_mp()) ) ;

    res.set_etat_qcq() ;

    res.set() = t1() - x ;	// Cmp - double

    return res ;

}


Qtenseur operator-(double x, const Qtenseur & t2) {
    
    return - (t2 - x) ; 
    
}


Qtenseur operator-(const Qtenseur & t1, int m) {

    return t1 - double(m) ; 
    
}


Qtenseur operator-(int m, const Qtenseur & t2) {

    return - (t2 - double(m)) ;     

}



			//****************//
			// MULTIPLICATION //
			//****************//



Qtenseur operator*(double x, const Qtenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    if ( (t.get_etat() == ETATZERO) || (x == double(1)) )
	return t ;
    else {
	Qtenseur res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
		    t.get_triad()) ;

	if ( x == double(0) )
	    res.set_etat_zero() ;
	else {
	    res.set_etat_qcq() ;
	    for (int i=0 ; i<res.get_n_comp() ; i++) {
		Itbl indices (res.donne_indices(i)) ;
		res.set(indices) = x*t(indices) ;
		}
	    }
	    return res ; 
	}
}


Qtenseur operator* (const Qtenseur& t, double x) {
    return x * t ;
}

Qtenseur operator*(int m, const Qtenseur& t) {
    return double(m) * t ; 
}


Qtenseur operator* (const Qtenseur& t, int m) {
    return double(m) * t ;
}


			//**********//
			// DIVISION //
			//**********//

Qtenseur operator/ (const Qtenseur& t1, const Qtenseur& t2) {
    
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t2.get_valence() == 0) ; // t2 doit etre un scalaire !
    assert(t1.get_mp() == t2.get_mp()) ;

    
    // Cas particuliers
    if (t2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Qtenseur / Qtenseur !" << endl ;
	abort() ; 
    }
    if (t1.get_etat() == ETATZERO) {
    	return t1 ;
    }

    // Cas general
    
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    Qtenseur res(*(t1.get_mp()), t1.get_valence(), t1.get_type_indice(), 
		t1.get_triad()) ;

    res.set_etat_qcq() ;
    for (int i=0 ; i<res.get_n_comp() ; i++) {
	Itbl indices (res.donne_indices(i)) ;
	res.set(indices) = t1(indices) / t2() ;	    // Cmp / Cmp
    }
    return res ;

}


Qtenseur operator/ (const Qtenseur& t, double x) {

    assert (t.get_etat() != ETATNONDEF) ;
 
    if ( x == double(0) ) {
	cout << "Division by 0 in Qtenseur / double !" << endl ;
	abort() ;
    }

    if ( (t.get_etat() == ETATZERO) || (x == double(1)) )
	return t ;
    else {
	Qtenseur res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
		    t.get_triad()) ;

	res.set_etat_qcq() ;
	for (int i=0 ; i<res.get_n_comp() ; i++) {
	    Itbl indices (t.donne_indices(i)) ;
	    res.set(indices) = t(indices) / x ;	    // Cmp / double
	}
	return res ; 
    }

}




Qtenseur operator/ (double x, const Qtenseur& t) {
    
    if (t.get_etat() == ETATZERO) {
	cout << "Division by 0 in double / Qtenseur !" << endl ;
	abort() ; 
    }
    
    assert (t.get_etat() == ETATQCQ) ;
    assert(t.get_valence() == 0) ;	// Utilisable que sur scalaire !
    
    Qtenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    res.set() = x / t() ;	// double / Cmp
    return res ;
}


Qtenseur operator/ (const Qtenseur& t, int m) {

    return t / double(m) ; 
}


Qtenseur operator/ (int m, const Qtenseur& t) {
    
    return double(m) / t ; 
}





