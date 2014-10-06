/*
 *  Arithmetical operations for class Scalar
 *
 */

/*
 *   Copyright (c) 2003-2005 Eric Gourgoulhon & Jerome Novak
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
 * Revision 1.9  2014/10/06 15:16:15  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2005/11/17 15:30:11  e_gourgoulhon
 * Added arithmetics with Mtbl.
 *
 * Revision 1.7  2004/07/06 13:36:29  j_novak
 * Added methods for desaliased product (operator |) only in r direction.
 *
 * Revision 1.6  2004/02/19 10:55:11  e_gourgoulhon
 * Treatment of case ETATUN in  double*Scalar, double/Scalar and
 * Scalar/double: added the copy of the spectral bases from the
 * input to the result.
 *
 * Revision 1.5  2003/11/03 22:35:45  e_gourgoulhon
 * Changed output comment when dzpuis conflict.
 *
 * Revision 1.4  2003/10/28 21:34:47  e_gourgoulhon
 * Corrected treatment of the case ETATUN in operator+=, operator-= and
 * operator*=.
 *
 * Revision 1.3  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.2  2003/10/01 13:04:44  e_gourgoulhon
 * The method Tensor::get_mp() returns now a reference (and not
 * a pointer) onto a mapping.
 *
 * Revision 1.1  2003/09/24 15:21:45  j_novak
 * New version
 *
 *
 * $Header$
 *
 */

// headers C
#include <cassert>
#include <cstdlib>

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
    assert( (ci.get_etat() == ETATQCQ) || (ci.get_etat() == ETATUN)) ;	
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
    if (c1.get_etat() == ETATUN) {
	return (c2 + double(1)) ;
    }
    if (c2.get_etat() == ETATUN) {
	return (c1 + double(1)) ;
    }
    assert(c1.get_etat() == ETATQCQ) ;
    assert(c2.get_etat() == ETATQCQ) ;
  
    // Cas general

    if ( c1.dz_nonzero() && c2.dz_nonzero() ) {
	if ( c1.get_dzpuis() != c2.get_dzpuis() ) {
	    cout << "Operation Scalar + Scalar: dzpuis conflict in the external " << endl;
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

// Scalar + Mtbl
// -------------
Scalar operator+(const Scalar& c1, const Mtbl& mi) {
    
    if ((c1.get_etat() == ETATNONDEF) || (mi.get_etat() == ETATNONDEF)) {
        cerr << "Undifined state in Scalar + Mtbl !" << endl ;
        abort() ; 
    }

    // Cas particuliers

    if (mi.get_etat() == ETATZERO) {
	return c1 ;
    }

    assert( c1.check_dzpuis(0) ) ; 
    
    Scalar resu(c1) ;
    
    if (c1.get_etat() == ETATZERO) {
	resu = mi ; 
    }
    else {
      if (c1.get_etat() == ETATUN) {
	resu = double(1) + mi ; 
      }
      else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.va = resu.va + mi ;
      }
    }
     
    resu.set_dzpuis(0) ; 
       
    return resu ;
}

// Mtbl + Scalar
// -------------
Scalar operator+(const Mtbl& mi, const Scalar& c1) {

    return c1 + mi ; 
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
    else {
      if (t1.get_etat() == ETATUN) {
	resu = double(1) + x ; 
      }
      else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.va = resu.va + x ;
      }
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
    if (c1.get_etat() == ETATUN) {
	return -(c2 - double(1)) ;
    }
    if (c2.get_etat() == ETATUN) {
	return (c1 - double(1)) ;
    }
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...

   // Cas general
   if ( c1.dz_nonzero() && c2.dz_nonzero() ) {
	if ( c1.get_dzpuis() != c2.get_dzpuis() ) {
	    cout << "Operation Scalar - Scalar : dzpuis conflict in the external " << endl;
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

// Scalar - Mtbl
// -------------
Scalar operator-(const Scalar& t1, const Mtbl& mi) {

    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    
    // Cas particuliers
    if (mi.get_etat() == ETATZERO) {
	return t1 ;
    }

    assert( t1.check_dzpuis(0) ) ; 
    
    Scalar resu(t1) ;
    
    if (t1.get_etat() == ETATZERO) {
      resu = - mi ; 
    }
    else{
      if (t1.get_etat() == ETATUN) {
	resu = double(1) - mi ; 
      }
      else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.va = resu.va - mi ;
      }
    } 
    resu.set_dzpuis(0) ; 

    return resu ;
}

// Mtbl - Scalar  
// -------------
Scalar operator-(const Mtbl& mi, const Scalar& t1) {

    return - (t1 - mi) ; 
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
      if (t1.get_etat() == ETATUN) {
	resu = double(1) - x ; 
      }
      else{ 
	assert(resu.get_etat() == ETATQCQ) ; // sinon ...
	resu.va = resu.va - x ;
      }
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
    if (c1.get_etat() == ETATUN)
	return c2 ;
    
    if (c2.get_etat() == ETATUN)
	return c1 ;
    
    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...
    
    // Protection
    assert( c1.get_mp() == c2.get_mp() ) ;
    
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
    if (c1.get_etat() == ETATUN)
	return c2 ;
    if (c2.get_etat() == ETATUN)
	return c1 ;

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

// Scalar | Scalar (multiplication with desaliasing in r)
// ------------------------------------------------------
Scalar operator|(const Scalar& c1, const Scalar& c2) {
    
    
    // Cas particuliers
    if ((c1.get_etat() == ETATZERO) || (c1.get_etat() == ETATNONDEF)){
	return c1 ;
    }
    if ((c2.get_etat() == ETATZERO)|| (c2.get_etat() == ETATNONDEF)) {
	return c2 ;
    }
    if (c1.get_etat() == ETATUN)
	return c2 ;
    if (c2.get_etat() == ETATUN)
	return c1 ;

    assert(c1.get_etat() == ETATQCQ) ;  // sinon...
    assert(c2.get_etat() == ETATQCQ) ;  // sinon...
    
    // Protection
    assert(c1.get_mp() == c2.get_mp()) ;
    
    // Cas general
    Scalar r( c1.get_mp() ) ;	    // Le resultat
    r.set_etat_qcq() ;
    r.va = c1.va | c2.va ;

    r.set_dzpuis( c1.get_dzpuis() + c2.get_dzpuis() ) ;
    
    // Termine
    return r ;
}


// Mtbl * Scalar
// -------------

Scalar operator*(const Mtbl& mi, const Scalar& c1) {

    // Particular cases 
    if ((c1.get_etat() == ETATZERO) || (c1.get_etat() == ETATNONDEF)) {
	return c1 ;
    }

    Scalar r(c1.get_mp()) ;
    if (c1.get_etat() == ETATUN) {
      r = mi ;
    }
    else {
      assert(c1.get_etat() == ETATQCQ) ;  // sinon...
      
      // Cas general
      r.set_dzpuis( c1.get_dzpuis() ) ;
      
      if ( mi.get_etat() == ETATZERO) {
	r.set_etat_zero() ;
      }
      else {
	r.set_etat_qcq() ;
	r.va = mi * c1.va ;
      }
    }

    // Termine
    return r ;
   
}

// Scalar * Mtbl
// -------------

Scalar operator*(const Scalar& c1, const Mtbl& mi) {

    return mi * c1 ; 
}

// double * Scalar
// ------------
Scalar operator*(double a, const Scalar& c1) {
    
    // Cas particuliers
    if ((c1.get_etat() == ETATZERO) || (c1.get_etat() == ETATNONDEF)) {
	return c1 ;
    }

    if (a == double(1))
      return c1 ;

    Scalar r(c1.get_mp()) ;
    if (c1.get_etat() == ETATUN) {
      r = a ;
      r.set_spectral_base(c1.get_spectral_va().get_base()) ; 
    }
    else {
      assert(c1.get_etat() == ETATQCQ) ;  // sinon...
      
      // Cas general
      r.set_dzpuis( c1.get_dzpuis() ) ;
      
      if ( a == double(0) ) {
	r.set_etat_zero() ;
      }
      else {
	r.set_etat_qcq() ;
	r.va = a * c1.va ;
      }
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
    if (c1.get_etat() == ETATUN)
	return double(1)/c2 ;
    if (c2.get_etat() == ETATUN)
	return c1 ;

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

// Scalar / Mtbl
// -------------
Scalar operator/(const Scalar& c1, const Mtbl& mi) {
  
    if (c1.get_etat() == ETATNONDEF) return c1 ;
	
    // Cas particuliers
    if ( mi.get_etat() == ETATZERO ) {
	cout << "Division by 0 in Scalar / Mtbl !" << endl ;
	abort() ;
    }
    if (c1.get_etat() == ETATZERO) {
	return c1 ;
    }
    Scalar r(c1.get_mp()) ;     // Le resultat

    if (c1.get_etat() == ETATUN) {
	r = double(1) / mi ;
    }
    else {
      assert(c1.get_etat() == ETATQCQ) ;  // sinon...

      r.set_etat_qcq() ;
      r.va = c1.va / mi ;

      r.set_dzpuis( c1.get_dzpuis() ) ;
    }
    // Termine
    return r ;
}


// Mtbl / Scalar
// -------------
Scalar operator/(const Mtbl& mi, const Scalar& c2) {
    
    if (c2.get_etat() == ETATNONDEF) 
	return c2 ;
	
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Mtbl / Scalar !" << endl ;
	abort() ; 
    }
    Scalar r(c2.get_mp()) ;     // Le resultat
    if (c2.get_etat() == ETATUN) {
      r = mi ;
    }
    else {
      assert(c2.get_etat() == ETATQCQ) ;  // sinon...
      
      r.set_dzpuis( - c2.get_dzpuis() ) ;
      
      if ( mi.get_etat() == ETATZERO ) {
	r.set_etat_zero() ;
      }
      else {
	r.set_etat_qcq() ;
	r.va = mi / c2.va ;
      }
    }

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
    Scalar r(c1.get_mp()) ;     // Le resultat

    if (c1.get_etat() == ETATUN) {
	r = double(1)/x ;
        r.set_spectral_base(c1.get_spectral_va().get_base()) ; 
    }
    else {
      assert(c1.get_etat() == ETATQCQ) ;  // sinon...

      r.set_etat_qcq() ;
      r.va = c1.va / x ;

      r.set_dzpuis( c1.get_dzpuis() ) ;
    }
    // Termine
    return r ;
}


// double / Scalar
// ------------
Scalar operator/(double x, const Scalar& c2) {
    
    if (c2.get_etat() == ETATNONDEF) 
	return c2 ;
	
    if (c2.get_etat() == ETATZERO) {
	cout << "Division by 0 in double / Scalar !" << endl ;
	abort() ; 
    }
    Scalar r(c2.get_mp()) ;     // Le resultat
    if (c2.get_etat() == ETATUN) {
      r = x ;
      r.set_spectral_base(c2.get_spectral_va().get_base()) ; 
    }
    else {
      assert(c2.get_etat() == ETATQCQ) ;  // sinon...
      
      r.set_dzpuis( - c2.get_dzpuis() ) ;
      
      if ( x == double(0) ) {
	r.set_etat_zero() ;
      }
      else {
	r.set_etat_qcq() ;
	r.va = x / c2.va ;
      }
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
    assert(mp == &(ci.get_mp()) ) ;	    // meme mapping
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
	    	cout << "Operation += Scalar : dzpuis conflict in the external " << endl;
	    	cout << " compactified domain ! " << endl ; 
	    	abort() ;
		}
    }
    
    if (etat == ETATZERO) {
		(*this) = ci ;
    }
    else {
		va += ci.va ;
    	if (etat == ETATUN) {
			etat = ETATQCQ ; 	// since the case ci.etat == ETATZERO
		}						//  has been treated above						
		
		assert(etat == ETATQCQ) ; 
	
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
    assert(mp == &(ci.get_mp()) ) ;	    // meme mapping
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
	    cout << "Operation -= Scalar : dzpuis conflict in the external " << endl;
	    cout << " compactified domain ! " << endl ; 
	    abort() ;
	}
    }
    

    if (etat == ETATZERO) {
		(*this) = -ci ;
    }
    else {
		va -= ci.va ;

    	if (etat == ETATUN) {
			etat = ETATQCQ ; 	// since the case ci.etat == ETATZERO
		}						//  has been treated above						
		
		assert(etat == ETATQCQ) ; 

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
    assert(mp == &(ci.get_mp()) ) ;	    // meme mapping
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
    
    if (ci.get_etat() == ETATUN) {
	return ;
    }
    
    if (etat == ETATUN) {
		operator=(ci) ; 
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
