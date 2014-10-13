/*
 * Arithmetical operations for class Itbl
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char itbl_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:01  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2002/10/16 14:36:38  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/11/17  16:04:54  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
// headers Lorene
#include "itbl.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

// + Itbl
// -----
namespace Lorene {
Itbl operator+(const Itbl& t1)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    return t1 ;
}

// - Itbl
// -----
Itbl operator-(const Itbl& t1)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    
    // Cas general
    Itbl r(t1.dim) ;		    // Itbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] = - (t1.t)[i] ;
    }
    return r ;
}

			//**********//
			// ADDITION //
			//**********//

// Itbl + Itbl
// ---------
Itbl operator+(const Itbl& t1, const Itbl& t2)
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_ndim() == t2.get_ndim()) ;
    for (int i=0 ; i<t1.get_ndim() ; i++) {
	assert( t1.get_dim(i) == t2.get_dim(i) ) ;
    }

    // Traitement des cas particuliers
    if (t1.get_etat() == ETATZERO) {
	return t2 ;
    }
    if (t2.get_etat() == ETATZERO) {
	return t1 ;
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    assert(t2.get_etat() == ETATQCQ) ;	// sinon...

    Itbl r(t1) ;	    // Itbl resultat
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] += (t2.t)[i] ;
    }
    
    // Termine
    return r ;
}

// Itbl + int
// ------------
Itbl operator+(const Itbl& t1, int x)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ( x == 0 ) {
	return t1 ;
    }
    
    // Cas general
    Itbl r(t1) ;		// Itbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] += x ;
    }
    return r ;
}

// int + Itbl
// ------------
Itbl operator+(int x, const Itbl& t1)
{
    return t1 + x ;
}


			//**************//
			// SOUSTRACTION //
			//**************//

// Itbl - Itbl
// ---------
Itbl operator-(const Itbl& t1, const Itbl& t2)
{

    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_ndim() == t2.get_ndim()) ;
    for (int i=0 ; i<t1.get_ndim() ; i++) {
	assert( t1.get_dim(i) == t2.get_dim(i) ) ;
    }

    // Traitement des cas particuliers
    if (t1.get_etat() == ETATZERO) {
	return -t2 ;
    }
    if (t2.get_etat() == ETATZERO) {
	return t1 ;
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    assert(t2.get_etat() == ETATQCQ) ;	// sinon...

    Itbl r(t1) ;	    // Itbl resultat
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] -= (t2.t)[i] ;
    }
    
    // Termine
    return r ;
}


// Itbl - int
// ------------
Itbl operator-(const Itbl& t1, int x)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ( x == 0 ) {
	return t1 ;
    }
    
    // Cas general
    Itbl r(t1) ;		// Itbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] -= x ;
    }
    return r ;
}


// int - Itbl
// ------------
Itbl operator-(int x, const Itbl& t1)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ( x == 0 ) {
	return -t1 ;
    }
    
    // Cas general
    Itbl r(t1) ;		// Itbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] -= x ;
    }
    return -r ;
}

			//****************//
			// MULTIPLICATION //
			//****************//

// Itbl * Itbl
// ---------
Itbl operator*(const Itbl& t1, const Itbl& t2)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_ndim() == t2.get_ndim()) ;
    for (int i=0 ; i<t1.get_ndim() ; i++) {
	assert( t1.get_dim(i) == t2.get_dim(i) ) ;
    }

    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    if (t2.get_etat() == ETATZERO) {
	return t2 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    assert(t2.get_etat() == ETATQCQ) ;	// sinon...

    Itbl r(t1) ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] *= (t2.t)[i] ;
    }
    
    // Termine
    return r ;
}

// Itbl * int
// ------------
Itbl operator*(const Itbl& t1, int x)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ((t1.get_etat() == ETATZERO) || ( x == 1 )) {
	return t1 ;
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Itbl r(t1) ;		    // Itbl resultat

    if (x == 0) {
	r.set_etat_zero() ;
    }
    else {
	for (int i=0 ; i<r.get_taille() ; i++) {
	    (r.t)[i] *= x ;
	}
    }
    
    // Termine
    return r ;
}

// int * Itbl
// ------------
Itbl operator*(int x, const Itbl& t1)
{
    return t1 * x ;
}


			//*******************//
			// operateurs +=,... //
			//*******************//

void Itbl::operator+=(const Itbl & ti) {

    // Protection
    assert(dim == ti.dim) ;
    assert(etat != ETATNONDEF) ;
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (ti.get_etat() == ETATZERO) {
    	return ;
    }
    
    // Cas general
    int n = get_taille() ;
    switch(etat) {
	case ETATZERO:
	set_etat_qcq() ;
	for (int i=0 ; i<n ; i++) {
    	    t[i] = ti.t[i] ;
	}
	break ;
	
	case ETATQCQ:
	for (int i=0 ; i<n ; i++) {
    	    t[i] += ti.t[i] ;
	}
	break ;
	
	default:
	cout << "etat inconnu " << __FILE__ << endl ;
	abort() ;
	break ;
    }
    
    // Termine
}

void Itbl::operator+=(int x) {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Cas particulier
    if ( x == 0 ) {
    	return ;
    }
    
    // Cas general
    int n = get_taille() ;
    switch(etat) {
	case ETATZERO:
	set_etat_qcq() ;
	for (int i=0 ; i<n ; i++) {
    	    t[i] = x ;
	}
	break ;
	
	case ETATQCQ:
	for (int i=0 ; i<n ; i++) {
    	    t[i] += x ;
	}
	break ;
	
	default:
	cout << "etat inconnu " << __FILE__ << endl ;
	abort() ;
	break ;
    }
    
    // Termine
}

void Itbl::operator-=(const Itbl & ti) {

    // Protection
    assert(dim == ti.dim) ;
    assert(etat != ETATNONDEF) ;
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (ti.get_etat() == ETATZERO) {
    	return ;
    }
    
    // Cas general
    int n = get_taille() ;
    switch(etat) {
	case ETATZERO:
	set_etat_qcq() ;
	for (int i=0 ; i<n ; i++) {
    	    t[i] = - ti.t[i] ;
	}
	break ;
	
	case ETATQCQ:
	for (int i=0 ; i<n ; i++) {
    	    t[i] -= ti.t[i] ;
	}
	break ;
	
	default:
	cout << "etat inconnu " << __FILE__ << endl ;
	abort() ;
	break ;
    }
    
    // Termine
}

void Itbl::operator-=(int x) {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Cas particulier
    if ( x == 0 ) {
    	return ;
    }
    
    // Cas general
    int n = get_taille() ;
    switch(etat) {
	case ETATZERO:
	set_etat_qcq() ;
	for (int i=0 ; i<n ; i++) {
    	    t[i] = - x ;
	}
	break ;
	
	case ETATQCQ:
	for (int i=0 ; i<n ; i++) {
    	    t[i] -= x ;
	}
	break ;
	
	default:
	cout << "etat inconnu " << __FILE__ << endl ;
	abort() ;
	break ;
    }
    
    // Termine
}

void Itbl::operator*=(const Itbl & ti) {

    // Protection
    assert(dim == ti.dim) ;
    assert(etat != ETATNONDEF) ;
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (etat == ETATZERO) {
    	return ;
    }
    if (ti.get_etat() == ETATZERO) {
    	set_etat_zero() ;
    	return ;
    }
    
    // Cas general
    assert(etat == ETATQCQ) ;
    for (int i=0 ; i<get_taille() ; i++) {
    	t[i] *= ti.t[i] ;
    }
    
    // Termine
}

void Itbl::operator*=(int x) {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Cas particulier
    if ( x == 0 ) {
    	set_etat_zero() ;
    	return ;
    }
    if (etat == ETATZERO) {
    	return ;
    }
    
    // Cas general
    assert(etat == ETATQCQ) ;
    for (int i=0 ; i<get_taille() ; i++) {
    	t[i] *= x ;
    }
    
    // Termine
}
}
