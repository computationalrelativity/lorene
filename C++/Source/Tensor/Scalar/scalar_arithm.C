/*
 *  Arithmetical operations for class Scalar
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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


char scalar_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/24 15:21:45  j_novak
 * New version
 *
 *
 * $Header$
 *
 */

// headers C
#include <assert.h>
#include <stdlib.h>

// headers Lorene
#include "tensor.h"
#include "type_parite.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

Scalar operator+(const Scalar & ci) {
    return ci ;
}

Scalar operator-(const Scalar & ci) {

    // Cas particulier
    if ((ci.get_etat() == ETATZERO) || (ci.get_etat() == ETATNONDEF)) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...
    Scalar r(ci.get_mp()) ;	// Scalar resultat
    r.set_etat_qcq() ;
    r.va = - ci.va ;    
    r.set_dzpuis( ci.get_dzpuis() ) ;
    
    // Termine
    return r ;
}

			//**********//
			// ADDITION //
			//**********//
// Scalar + Scalar
// ---------
Scalar operator+(const Scalar & c1, const Scalar & c2) {
    
    if (c1.get_etat() == ETATNONDEF) 
	return c1 ;
    if (c2.get_etat() == ETATNONDEF) 
	return c2 ;
    assert(c1.get_mp() == c2.get_mp()) ;
    
    // Cas particuliers
    if (c1.get_etat() == ETATZERO) {
	return c2 ;
    }
    if (c2.get_etat() == ETATZERO) {
	return c1 ;
    }
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...
  
    // Cas general

    if ( c1.dz_nonzero() && c2.dz_nonzero() ) {
	if ( c1.get_dzpuis() != c2.get_dzpuis() ) {
	    cout << "Operation Scalar + Scalar forbidden in the external " << endl;
	    cout << " compactified domain ! " << endl ; 
	    abort() ;
	}
    }
    
    Scalar r(c1) ;	    // Le resultat
    r.va += c2.va ;
    
    if (c1.dz_nonzero()) {
	r.set_dzpuis( c1.get_dzpuis() ) ; 
    }
    else{
	r.set_dzpuis( c2.get_dzpuis() ) ; 	
    }

    // Termine
    return r ;
}

// Scalar + double
// ------------
Scalar operator+(const Scalar& t1, double x)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (x == double(0)) {
	return t1 ;
    }

    assert( t1.check_dzpuis(0) ) ; 
     
    Scalar resu(t1) ;
    
    if (t1.get_etat() == ETATZERO) {
	resu = x ; 
    }
    else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.va = resu.va + x ;
    }
     
    resu.set_dzpuis(0) ; 
       
    return resu ;
}

// double + Scalar
// ------------
Scalar operator+(double x, const Scalar& t1)	   
{
    return t1 + x ;
}

// Scalar + int
// ---------
Scalar operator+(const Scalar& t1, int m)	    
{
    return t1 + double(m) ;
}

// int + Scalar
// ---------
Scalar operator+(int m, const Scalar& t1)	   
{
    return t1 + double(m) ;
}





			//**************//
			// SOUSTRACTION //
			//**************//

// Scalar - Scalar
// ---------
Scalar operator-(const Scalar & c1, const Scalar & c2) {
    
    if (c1.get_etat() == ETATNONDEF) 
	return c1 ;
    if (c2.get_etat() == ETATNONDEF) 
	return c2 ;
    
    assert(c1.get_mp() == c2.get_mp()) ;
    
    // Cas particuliers
    if (c1.get_etat() == ETATZERO) {
	return -c2 ;
    }
    if (c2.get_etat() == ETATZERO) {
	return c1 ;
    }
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...

   // Cas general
   if ( c1.dz_nonzero() && c2.dz_nonzero() ) {
	if ( c1.get_dzpuis() != c2.get_dzpuis() ) {
	    cout << "Operation Scalar - Scalar forbidden in the external " << endl;
	    cout << " compactified domain ! " << endl ; 
	    abort() ;
	}
    }
    
    Scalar r(c1) ;	    // Le resultat
    r.va -= c2.va ;
    
    if (c1.dz_nonzero()) {
	r.set_dzpuis( c1.get_dzpuis() ) ; 
    }
    else{
	r.set_dzpuis( c2.get_dzpuis() ) ; 	
    }
        
    // Termine
    return r ;
}

// Scalar - double
// ------------
Scalar operator-(const Scalar& t1, double x)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (x == double(0)) {
	return t1 ;
    }

    assert( t1.check_dzpuis(0) ) ; 

    Scalar resu(t1) ;
    
    if (t1.get_etat() == ETATZERO) {
	resu = - x ; 
    }
    else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.va = resu.va - x ;
    }
        
    resu.set_dzpuis(0) ; 

    return resu ;
}

// double - Scalar
// ------------
Scalar operator-(double x, const Scalar& t1)	    
{
    return - (t1 - x) ;
}

// Scalar - int
// ---------
Scalar operator-(const Scalar& t1, int m)	    
{
    return t1 - double(m) ;
}

// int - Scalar
// ---------
Scalar operator-(int m, const Scalar& t1)	    
{
    return double(m) - t1 ;
}






			//****************//
			// MULTIPLICATION //
			//****************//

// Scalar * Scalar
// ---------
Scalar operator*(const Scalar& c1, const Scalar& c2) {
    
    
    // Cas particuliers
    if ((c1.get_etat() == ETATZERO) || (c1.get_etat() == ETATNONDEF)){
	return c1 ;
    }
    if ((c2.get_etat() == ETATZERO)|| (c2.get_etat() == ETATNONDEF)) {
	return c2 ;
    }
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...
    
    // Protection
    assert(*(c1.get_mp()) == *(c2.get_mp())) ;
    
    // Cas general
    Scalar r(c1) ;	    // Le resultat
    r.va *= c2.va ;

    r.set_dzpuis( c1.get_dzpuis() + c2.get_dzpuis() ) ;
    
    // Termine
    return r ;
}

// Scalar % Scalar (multiplication with desaliasing)
// -------------------------------------------
Scalar operator%(const Scalar& c1, const Scalar& c2) {
    
    
    // Cas particuliers
    if ((c1.get_etat() == ETATZERO) || (c1.get_etat() == ETATNONDEF)){
	return c1 ;
    }
    if ((c2.get_etat() == ETATZERO)|| (c2.get_etat() == ETATNONDEF)) {
	return c2 ;
    }
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...
    
    // Protection
    assert(c1.get_mp() == c2.get_mp()) ;
    
    // Cas general
    Scalar r( c1.get_mp() ) ;	    // Le resultat
    r.set_etat_qcq() ;
    r.va = c1.va % c2.va ;

    r.set_dzpuis( c1.get_dzpuis() + c2.get_dzpuis() ) ;
    
    // Termine
    return r ;
}




// double * Scalar
// ------------
Scalar operator*(double a, const Scalar& c1) {
    
    // Cas particuliers
    if ((c1.get_etat() == ETATZERO) || (c1.get_etat() == ETATNONDEF)) {
	return c1 ;
    }

    assert(c1.get_etat() == ETATQCQ) ;  // sinon...

    // Cas general
    Scalar r(c1.get_mp()) ;
    r.set_dzpuis( c1.get_dzpuis() ) ;
    
    if ( a == double(0) ) {
	r.set_etat_zero() ;
    }
    else {
	r.set_etat_qcq() ;
	r.va = a * c1.va ;
    }
    

    // Termine
    return r ;
}


// Scalar * double
// ------------
Scalar operator*(const Scalar& t1, double x)	    
{
    return x * t1 ;
}

// Scalar * int
// ---------
Scalar operator*(const Scalar& t1, int m)	    
{
    return t1 * double(m) ;
}

// int * Scalar
// ---------
Scalar operator*(int m, const Scalar& t1)	    
{
    return double(m) * t1 ;
}







			//**********//
			// DIVISION //
			//**********//


// Scalar / Scalar
// ---------
Scalar operator/(const Scalar& c1, const Scalar& c2) {
    
    // Protections
    assert(c1.get_etat() != ETATNONDEF) ;
    assert(c2.get_etat() != ETATNONDEF) ;
    assert(c1.get_mp() == c2.get_mp()) ;
    
    // Cas particuliers
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Scalar / Scalar !" << endl ;
	abort() ; 
    }
    if (c1.get_etat() == ETATZERO) {
    	return c1 ;
    }

    // Cas general
    
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...

    Scalar r(c1.get_mp()) ;	    // Le resultat

    r.set_etat_qcq() ;
    r.va = c1.va / c2.va ;
    
    r.set_dzpuis( c1.get_dzpuis() - c2.get_dzpuis() ) ;

    // Termine
    return r ;
}

// Scalar / double
// -------------
Scalar operator/(const Scalar& c1, double x) {
    
    if (c1.get_etat() == ETATNONDEF) 
	return c1 ;
	
    // Cas particuliers
    if ( x == double(0) ) {
	cout << "Division by 0 in Scalar / double !" << endl ;
	abort() ;
    }
    if (c1.get_etat() == ETATZERO) {
	return c1 ;
    }
    
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...

    Scalar r(c1.get_mp()) ;     // Le resultat
 
    r.set_etat_qcq() ;
    r.va = c1.va / x ;

    r.set_dzpuis( c1.get_dzpuis() ) ;

    // Termine
    return r ;
}


// double / Scalar
// ------------
Scalar operator/(double x, const Scalar& c2) {
    
    if (c2.get_etat() == ETATNONDEF) 
	return c2 ;
	
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Scalar / Scalar !" << endl ;
	abort() ; 
    }
    
   
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...

    Scalar r(c2.get_mp()) ;     // Le resultat
    r.set_dzpuis( - c2.get_dzpuis() ) ;
 
    if ( x == double(0) ) {
	r.set_etat_zero() ;
    }
    else {
	r.set_etat_qcq() ;
	r.va = x / c2.va ;
    }

    // Termine
    return r ;
}


// Scalar / int
// ---------
Scalar operator/(const Scalar& c1, int m) {
    
    return c1 / double(m) ; 

}


// int / Scalar
// ---------
Scalar operator/(int m, const Scalar& c2) {

    return double(m) / c2 ;

}

			//*******************//
			// operateurs +=,... //
			//*******************//

//---------
//  += Scalar
//---------

void Scalar::operator+=(const Scalar & ci) {
    
    // Protection
    assert(mp == ci.get_mp()) ;	    // meme mapping
    if (etat == ETATNONDEF) 
	return ;
 
    // Cas particulier
    if (ci.get_etat() == ETATZERO) {
	return ;
    }
    
    if (ci.get_etat() == ETATNONDEF) {
	set_etat_nondef() ;
	return ;
    }
        
    // Cas general
    

    if ( dz_nonzero() && ci.dz_nonzero() ) {
	if ( dzpuis != ci.dzpuis ) {
	    cout << "Operation += Scalar forbidden in the external " << endl;
	    cout << " compactified domain ! " << endl ; 
	    abort() ;
	}
    }
    
    if (etat == ETATZERO) {
	(*this) = ci ;
    }
    else {
	va += ci.va ;
    
	if( ci.dz_nonzero() ) {
	    set_dzpuis(ci.dzpuis) ; 
	}
    }
    // Menage (a ne faire qu'a la fin seulement)
    del_deriv() ;

    
}

//---------
//  -= Scalar
//---------

void Scalar::operator-=(const Scalar & ci) {
    
    // Protection
    assert(mp == ci.get_mp()) ;	    // meme mapping
    if (etat == ETATNONDEF) 
	return ;
 
    // Cas particulier
    if (ci.get_etat() == ETATZERO) {
	return ;
    }
    
    if (ci.get_etat() == ETATNONDEF) {
	set_etat_nondef() ;
	return ;
    }
    
    // Cas general
    if ( dz_nonzero() && ci.dz_nonzero() ) {
	if ( dzpuis != ci.dzpuis ) {
	    cout << "Operation -= Scalar forbidden in the external " << endl;
	    cout << " compactified domain ! " << endl ; 
	    abort() ;
	}
    }
    

    if (etat == ETATZERO) {
	(*this) = -ci ;
    }
    else {
	va -= ci.va ;

	if( ci.dz_nonzero() ) {
	    set_dzpuis(ci.dzpuis) ; 
	}
    }
    // Menage (a ne faire qu'a la fin seulement)
    del_deriv() ;
}

//---------
//  *= Scalar
//---------

void Scalar::operator*=(const Scalar & ci) {
    
    // Protection
    assert(mp == ci.get_mp()) ;	    // meme mapping
    if (etat == ETATNONDEF) 
	return ;
 
    // Cas particulier
    if (ci.get_etat() == ETATZERO) {
	set_etat_zero() ;
	return ;
    }
    
    if (etat == ETATZERO) {
	return ; 
    }
    
    if (ci.get_etat() == ETATNONDEF) {
	set_etat_nondef() ;
	return ;
    }
        
    // Cas general
    
    assert(etat == ETATQCQ) ; // sinon....
    
    va *= ci.va ;
 
    dzpuis += ci.dzpuis ;    

    // Menage (a ne faire qu'a la fin seulement)
    del_deriv() ;

}
