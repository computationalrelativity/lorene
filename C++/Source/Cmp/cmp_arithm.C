/*
 *  Arithmetical operations for class Cmp
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char cmp_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:13:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2003/06/20 13:59:59  f_limousin
 * L'assert pour le mapping est realise a partir du mapping lui meme et non a partir du pointeur sur le mapping.
 *
 * Revision 1.2  2002/10/16 14:36:34  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.11  2001/05/26  15:07:47  eric
 * Ajout de operator% : multiplication de deux Cmp avec desaliasage
 *
 * Revision 2.10  2000/02/15  13:39:57  phil
 * correction -= et +=
 *
 * Revision 2.9  2000/01/28  16:34:57  eric
 * Utilisation des nouvelles fonctions Cmp::check_dzpuis et Cmp::dz_nonzero
 * dans les tests sur dzpuis.
 *
 * Revision 2.8  1999/12/10  16:32:52  eric
 * Dans l'arithmetique membre (+=, -=, *=), on n'appelle desormais
 *  del_deriv() que tout a la fin.
 *
 * Revision 2.7  1999/11/26  14:23:38  eric
 * Ajout du membre dzpuis.
 *
 * Revision 2.6  1999/11/12  17:08:35  eric
 * Ajout de la partie manquante de l'arithmetique.
 *
 * Revision 2.5  1999/10/28  13:23:43  phil
 * Deverouillage des ETATNONDEF
 *
 * Revision 2.4  1999/10/27  08:45:21  eric
 * ntroduction du membre Valeur va.
 *
 * Revision 2.3  1999/10/22  08:14:58  eric
 * La fonction annule() est rebaptisee annule_hard().
 *
 * Revision 2.2  1999/09/14  17:19:44  phil
 * ajout de Cmp operator* (double, const Cmp&)
 *
 * Revision 2.1  1999/03/03  11:13:56  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers C
#include <cassert>
#include <cstdlib>

// headers Lorene
#include "cmp.h"
#include "type_parite.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

Cmp operator+(const Cmp & ci) {
    return ci ;
}

Cmp operator-(const Cmp & ci) {

    // Cas particulier
    if ((ci.get_etat() == ETATZERO) || (ci.get_etat() == ETATNONDEF)) {
	return ci ;
    }
    
    // Cas general
    assert(ci.get_etat() == ETATQCQ) ;	// sinon...
    Cmp r(ci.get_mp()) ;	// Cmp resultat
    r.set_etat_qcq() ;
    r.va = - ci.va ;    
    r.set_dzpuis( ci.get_dzpuis() ) ;
    
    // Termine
    return r ;
}

			//**********//
			// ADDITION //
			//**********//
// Cmp + Cmp
// ---------
Cmp operator+(const Cmp & c1, const Cmp & c2) {
    
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
	    cout << "Operation Cmp + Cmp forbidden in the external " << endl;
	    cout << " compactified domain ! " << endl ; 
	    abort() ;
	}
    }
    
    Cmp r(c1) ;	    // Le resultat
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

// Cmp + double
// ------------
Cmp operator+(const Cmp& t1, double x)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (x == double(0)) {
	return t1 ;
    }

    assert( t1.check_dzpuis(0) ) ; 
     
    Cmp resu(t1) ;
    
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

// double + Cmp
// ------------
Cmp operator+(double x, const Cmp& t1)	   
{
    return t1 + x ;
}

// Cmp + int
// ---------
Cmp operator+(const Cmp& t1, int m)	    
{
    return t1 + double(m) ;
}

// int + Cmp
// ---------
Cmp operator+(int m, const Cmp& t1)	   
{
    return t1 + double(m) ;
}





			//**************//
			// SOUSTRACTION //
			//**************//

// Cmp - Cmp
// ---------
Cmp operator-(const Cmp & c1, const Cmp & c2) {
    
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
	    cout << "Operation Cmp - Cmp forbidden in the external " << endl;
	    cout << " compactified domain ! " << endl ; 
	    abort() ;
	}
    }
    
    Cmp r(c1) ;	    // Le resultat
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

// Cmp - double
// ------------
Cmp operator-(const Cmp& t1, double x)	   
{
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (x == double(0)) {
	return t1 ;
    }

    assert( t1.check_dzpuis(0) ) ; 

    Cmp resu(t1) ;
    
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

// double - Cmp
// ------------
Cmp operator-(double x, const Cmp& t1)	    
{
    return - (t1 - x) ;
}

// Cmp - int
// ---------
Cmp operator-(const Cmp& t1, int m)	    
{
    return t1 - double(m) ;
}

// int - Cmp
// ---------
Cmp operator-(int m, const Cmp& t1)	    
{
    return double(m) - t1 ;
}






			//****************//
			// MULTIPLICATION //
			//****************//

// Cmp * Cmp
// ---------
Cmp operator*(const Cmp& c1, const Cmp& c2) {
    
    
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
    Cmp r(c1) ;	    // Le resultat
    r.va *= c2.va ;

    r.set_dzpuis( c1.get_dzpuis() + c2.get_dzpuis() ) ;
    
    // Termine
    return r ;
}

// Cmp % Cmp (multiplication with desaliasing)
// -------------------------------------------
Cmp operator%(const Cmp& c1, const Cmp& c2) {
    
    
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
    Cmp r( c1.get_mp() ) ;	    // Le resultat
    r.set_etat_qcq() ;
    r.va = c1.va % c2.va ;

    r.set_dzpuis( c1.get_dzpuis() + c2.get_dzpuis() ) ;
    
    // Termine
    return r ;
}




// double * Cmp
// ------------
Cmp operator*(double a, const Cmp& c1) {
    
    // Cas particuliers
    if ((c1.get_etat() == ETATZERO) || (c1.get_etat() == ETATNONDEF)) {
	return c1 ;
    }

    assert(c1.get_etat() == ETATQCQ) ;  // sinon...

    // Cas general
    Cmp r(c1.get_mp()) ;
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


// Cmp * double
// ------------
Cmp operator*(const Cmp& t1, double x)	    
{
    return x * t1 ;
}

// Cmp * int
// ---------
Cmp operator*(const Cmp& t1, int m)	    
{
    return t1 * double(m) ;
}

// int * Cmp
// ---------
Cmp operator*(int m, const Cmp& t1)	    
{
    return double(m) * t1 ;
}







			//**********//
			// DIVISION //
			//**********//


// Cmp / Cmp
// ---------
Cmp operator/(const Cmp& c1, const Cmp& c2) {
    
    // Protections
    assert(c1.get_etat() != ETATNONDEF) ;
    assert(c2.get_etat() != ETATNONDEF) ;
    assert(c1.get_mp() == c2.get_mp()) ;
    
    // Cas particuliers
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Cmp / Cmp !" << endl ;
	abort() ; 
    }
    if (c1.get_etat() == ETATZERO) {
    	return c1 ;
    }

    // Cas general
    
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...

    Cmp r(c1.get_mp()) ;	    // Le resultat

    r.set_etat_qcq() ;
    r.va = c1.va / c2.va ;
    
    r.set_dzpuis( c1.get_dzpuis() - c2.get_dzpuis() ) ;

    // Termine
    return r ;
}

// Cmp / double
// -------------
Cmp operator/(const Cmp& c1, double x) {
    
    if (c1.get_etat() == ETATNONDEF) 
	return c1 ;
	
    // Cas particuliers
    if ( x == double(0) ) {
	cout << "Division by 0 in Cmp / double !" << endl ;
	abort() ;
    }
    if (c1.get_etat() == ETATZERO) {
	return c1 ;
    }
    
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...

    Cmp r(c1.get_mp()) ;     // Le resultat
 
    r.set_etat_qcq() ;
    r.va = c1.va / x ;

    r.set_dzpuis( c1.get_dzpuis() ) ;

    // Termine
    return r ;
}


// double / Cmp
// ------------
Cmp operator/(double x, const Cmp& c2) {
    
    if (c2.get_etat() == ETATNONDEF) 
	return c2 ;
	
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Cmp / Cmp !" << endl ;
	abort() ; 
    }
    
   
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...

    Cmp r(c2.get_mp()) ;     // Le resultat
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


// Cmp / int
// ---------
Cmp operator/(const Cmp& c1, int m) {
    
    return c1 / double(m) ; 

}


// int / Cmp
// ---------
Cmp operator/(int m, const Cmp& c2) {

    return double(m) / c2 ;

}

			//*******************//
			// operateurs +=,... //
			//*******************//

//---------
//  += Cmp
//---------

void Cmp::operator+=(const Cmp & ci) {
    
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
	    cout << "Operation += Cmp forbidden in the external " << endl;
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
//  -= Cmp
//---------

void Cmp::operator-=(const Cmp & ci) {
    
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
	    cout << "Operation -= Cmp forbidden in the external " << endl;
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
//  *= Cmp
//---------

void Cmp::operator*=(const Cmp & ci) {
    
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
