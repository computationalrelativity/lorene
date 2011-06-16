/*
 * Arithmetical operations for class Tbl
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


char tbl_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2011/06/16 10:48:28  j_novak
 * Minor modif.
 *
 * Revision 1.2  2002/10/16 14:37:14  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/09/27  14:20:39  eric
 * Remplacement du text x == 0. par x == double(0)
 *  dans la multiplication par un double.
 *
 * Revision 2.3  1999/11/15  16:36:55  eric
 * Tbl::dim est desormais un Dim_tbl et non plus un Dim_tbl*.
 *
 * Revision 2.2  1999/10/01  10:09:34  eric
 * 0 -> double(0)
 *
 * Revision 2.1  1999/09/24  14:23:55  eric
 * Changement de prototypes.
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// headers Lorene
#include "tbl.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

// + Tbl
// -----
Tbl operator+(const Tbl& t1)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    return t1 ;
}

// - Tbl
// -----
Tbl operator-(const Tbl& t1)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    
    // Cas general
    Tbl r(t1.dim) ;		    // Tbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] = - (t1.t)[i] ;
    }
    return r ;
}

			//**********//
			// ADDITION //
			//**********//

// Tbl + Tbl
// ---------
Tbl operator+(const Tbl& t1, const Tbl& t2)
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

    Tbl r(t1) ;	    // Tbl resultat
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] += (t2.t)[i] ;
    }
    
    // Termine
    return r ;
}

// Tbl + double
// ------------
Tbl operator+(const Tbl& t1, double x)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ( x == double(0) ) {
	return t1 ;
    }
    
    // Cas general
    Tbl r(t1) ;		// Tbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] += x ;
    }
    return r ;
}

// double + Tbl
// ------------
Tbl operator+(double x, const Tbl& t1)
{
    return t1 + x ;
}

// Tbl + int
// ---------
Tbl operator+(const Tbl& t1, int n)
{
    return t1 + double(n) ;
}

// int + Tbl
// ---------
Tbl operator+(int n, const Tbl& t1)
{
    return t1 + double(n) ;
}


			//**************//
			// SOUSTRACTION //
			//**************//

// Tbl - Tbl
// ---------
Tbl operator-(const Tbl& t1, const Tbl& t2)
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

    Tbl r(t1) ;	    // Tbl resultat
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] -= (t2.t)[i] ;
    }
    
    // Termine
    return r ;
}


// Tbl - double
// ------------
Tbl operator-(const Tbl& t1, double x)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ( x == double(0) ) {
	return t1 ;
    }
    
    // Cas general
    Tbl r(t1) ;		// Tbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] -= x ;
    }
    return r ;
}

// Tbl - int
// ---------
Tbl operator-(const Tbl& t1, int n)
{
    return t1 - double(n) ;
}

// double - Tbl
// ------------
Tbl operator-(double x, const Tbl& t1)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ( x == double(0) ) {
	return -t1 ;
    }
    
    // Cas general
    Tbl r(t1) ;		// Tbl resultat
    r.set_etat_qcq() ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] -= x ;
    }
    return -r ;
}

// int - Tbl
// ---------
Tbl operator-(int n, const Tbl& t1)
{
    return double(n) - t1 ;
}

			//****************//
			// MULTIPLICATION //
			//****************//

// Tbl * Tbl
// ---------
Tbl operator*(const Tbl& t1, const Tbl& t2)
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

    Tbl r(t1) ;
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] *= (t2.t)[i] ;
    }
    
    // Termine
    return r ;
}

// Tbl * double
// ------------
Tbl operator*(const Tbl& t1, double x)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particulier
    if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
	return t1 ;
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Tbl r(t1) ;		    // Tbl resultat

    if (x == double(0)) {
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

// double * Tbl
// ------------
Tbl operator*(double x, const Tbl& t1)
{
    return t1 * x ;
}

// Tbl * int
// ---------
Tbl operator*(const Tbl& t1, int n)
{
    return t1 * double(n) ;
}

// int * Tbl
// ---------
Tbl operator*(int n, const Tbl& t1)
{
    return t1 * double(n) ;
}

			//**********//
			// DIVISION //
			//**********//

// Tbl / Tbl
// ---------
Tbl operator/(const Tbl& t1, const Tbl& t2)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t1.get_ndim() == t2.get_ndim()) ;
    for (int i=0 ; i<t1.get_ndim() ; i++) {
	assert( t1.get_dim(i) == t2.get_dim(i) ) ;
    }

    // Cas particuliers
    if (t2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Tbl/Tbl !" << endl ;
	abort() ; 
    }
    if (t1.get_etat() == ETATZERO) {
	return t1 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...
    assert(t2.get_etat() == ETATQCQ) ;	// sinon...

    Tbl r(t1) ;		    // Tbl resultat
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] /= (t2.t)[i] ;
    }
    
    // Termine
    return r ;
}

// Tbl / double
// ------------
Tbl operator/(const Tbl& t1, double x)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;
    if ( x == double(0) ) {
	cout << "Division by 0 in Tbl/double !" << endl ;
	abort() ;
    }
    
    // Cas particulier
    if ((t1.get_etat() == ETATZERO) || ( x == double(1) )) {
	return t1 ;
    }
    
    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Tbl r(t1) ;		    // Tbl resultat
    for (int i=0 ; i<r.get_taille() ; i++) {
	(r.t)[i] /= x ;
    }
    return r ;
}

// Tbl / int
// ---------
Tbl operator/(const Tbl& t1, int n)
{
    return t1 / double(n) ;
}

// double / Tbl
// ------------
Tbl operator/(double x, const Tbl& t1)
{
    // Protection
    assert(t1.get_etat() != ETATNONDEF) ;

    // Cas particuliers
    if (t1.get_etat() == ETATZERO) {
	cout << "Division by 0 in double/Tbl !" << endl ;
	abort() ; 
    }

    // Cas general
    assert(t1.get_etat() == ETATQCQ) ;	// sinon...

    Tbl r(t1.dim) ;		    	// Tbl resultat, a priori NONDEF
    
    if ( x == double(0) ) {
	r.set_etat_zero() ;
    }
    else {
    	r.set_etat_qcq() ;
	for (int i=0 ; i<r.get_taille() ; i++) {
	    (r.t)[i] = x / (t1.t)[i] ;
	}
    }
    
    // Termine
    return r ;
}

// int / Tbl
// ---------
Tbl operator/(int n, const Tbl& t1)
{
    return double(n) / t1 ;
}

			//*******************//
			// operateurs +=,... //
			//*******************//

void Tbl::operator+=(const Tbl & ti) {

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

void Tbl::operator+=(double x) {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Cas particulier
    if ( x == double(0) ) {
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

void Tbl::operator-=(const Tbl & ti) {

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

void Tbl::operator-=(double x) {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Cas particulier
    if ( x == double(0) ) {
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

void Tbl::operator*=(const Tbl & ti) {

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
    int n = get_taille() ;
    for (int i=0 ; i<n ; i++) {
    	t[i] *= ti.t[i] ;
    }
    
    // Termine
}

void Tbl::operator*=(double x) {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Cas particulier
    if ( x == double(0) ) {
    	set_etat_zero() ;
    	return ;
    }
    if (etat == ETATZERO) {
    	return ;
    }
    
    // Cas general
    int n = get_taille() ;
    assert(etat == ETATQCQ) ;
    for (int i=0 ; i<n ; i++) {
    	t[i] *= x ;
    }
    
    // Termine
}

void Tbl::operator/=(const Tbl & ti) {

    // Protection
    assert(dim == ti.dim) ;
    assert(etat != ETATNONDEF) ;
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (ti.get_etat() == ETATZERO) {
    	cout << "Division by 0 in Tbl::operator/=(const Tbl &) !" << endl ;
	abort() ;
    }
    if (etat == ETATZERO) {
    	return ;
    }
    
    // Cas general
    assert(etat == ETATQCQ) ;
    assert(ti.get_etat() == ETATQCQ) ;
    int n = get_taille() ;
    for (int i=0 ; i<n ; i++) {
    	t[i] /= ti.t[i] ;
    }
    
    // Termine
}

void Tbl::operator/=(double x) {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Cas particulier
    if ( x == double(0) ) {
    	cout << "Division by 0 in Tbl::operator/=(double ) !" << endl ;
	abort() ;
    }
    if (etat == ETATZERO) {
	return ;
    }
    
    // Cas general
    assert(etat == ETATQCQ) ;
    int n = get_taille() ;
    for (int i=0 ; i<n ; i++) {
    	t[i] /= x ;
    }
    
    // Termine
}

