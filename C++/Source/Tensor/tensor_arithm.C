/*
 *  Arithmetics functions for the Tensor class.
 *
 *  These functions are not member functions of the Tensor class.
 *
 *  (see file tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *   Copyright (c) 1999-2001 Philippe Grandclement (Tenseur version)
 *   Copyright (c) 2000-2001 Eric Gourgoulhon      (Tenseur version)
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:18:17  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:53:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:20  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2004/01/08 09:24:11  e_gourgoulhon
 * Added arithmetics common to Scalar and Tensor.
 * Corrected treatment ETATUN in Tensor / Scalar.
 *
 * Revision 1.2  2003/10/01 13:04:44  e_gourgoulhon
 * The method Tensor::get_mp() returns now a reference (and not
 * a pointer) onto a mapping.
 *
 * Revision 1.1  2003/09/26 14:33:53  j_novak
 * Arithmetic functions for the class Tensor
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>
#include <cmath>

// Headers Lorene
#include "tensor.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

namespace Lorene {
Tensor operator+(const Tensor & t) {

    return t ; 

}

Tensor operator-(const Tensor & t) {
    
  Tensor res(t.get_mp(), t.get_valence(), t.get_index_type(), 
		    t.get_triad()) ;

  for (int i=0 ; i<res.get_n_comp() ; i++) {
    Itbl ind (res.indices(i)) ;    
    res.set(ind) = -t(ind) ;
  }
  return res ;

}

			//**********//
			// ADDITION //
			//**********//

Tensor operator+(const Tensor & t1, const Tensor & t2) {
    
    assert (t1.get_valence() == t2.get_valence()) ;
    assert (t1.get_mp() == t2.get_mp()) ;
    if (t1.get_valence() != 0) {
	assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    for (int i=0 ; i<t1.get_valence() ; i++)
	assert(t1.get_index_type(i) == t2.get_index_type(i)) ;

    Tensor res(t1.get_mp(), t1.get_valence(), t1.get_index_type(), 
			t1.get_triad()) ;

    for (int i=0 ; i<res.get_n_comp() ; i++) {
      Itbl ind (res.indices(i)) ;
      res.set(ind) = t1(ind) + t2(ind) ;
    }
    return res ;

}


Scalar operator+(const Tensor& t1, const Scalar& t2) {
    
    assert (t1.get_valence() == 0) ;
    assert (t1.get_mp() == t2.get_mp()) ;

    return  *(t1.cmp[0]) + t2 ; 

}

Scalar operator+(const Scalar& t1, const Tensor& t2) {
    
    assert (t2.get_valence() == 0) ;
    assert (t1.get_mp() == t2.get_mp()) ;

    return  t1 + *(t2.cmp[0]) ; 

}

			//**************//
			// SOUSTRACTION //
			//**************//

Tensor operator-(const Tensor & t1, const Tensor & t2) {

    assert (t1.get_valence() == t2.get_valence()) ;
    assert (t1.get_mp() == t2.get_mp()) ;
    if (t1.get_valence() != 0) {
	assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    for (int i=0 ; i<t1.get_valence() ; i++)
	assert(t1.get_index_type(i) == t2.get_index_type(i)) ;

    Tensor res(t1.get_mp(), t1.get_valence(), t1.get_index_type(), 
			t1.get_triad()) ;

    for (int i=0 ; i<res.get_n_comp() ; i++) {
      Itbl ind (res.indices(i)) ;
      res.set(ind) = t1(ind) - t2(ind) ;
    }
    return res ;

}

Scalar operator-(const Tensor& t1, const Scalar& t2) {
    
    assert (t1.get_valence() == 0) ;
    assert (t1.get_mp() == t2.get_mp()) ;

    return  *(t1.cmp[0]) - t2 ; 

}

Scalar operator-(const Scalar& t1, const Tensor& t2) {
    
    assert (t2.get_valence() == 0) ;
    assert (t1.get_mp() == t2.get_mp()) ;

    return  t1 - *(t2.cmp[0]) ; 

}



			//****************//
			// MULTIPLICATION //
			//****************//

Tensor operator*(const Scalar& t1, const Tensor& t2) {
   
    assert (&(t1.get_mp()) == &(t2.get_mp())) ;
        
    if (t1.get_etat() == ETATUN) return t2 ;
    
    Tensor res(t2.get_mp(), t2.get_valence(), t2.get_index_type(), 
               t2.get_triad()) ;
    
    for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
        Itbl ind = res.indices(ic) ;
        res.set(ind) = t1 * t2(ind) ;
    }
    
    return res ;
}


Tensor operator*(const Tensor& t2, const Scalar& t1) {
       
    return t1*t2 ;
}



Tensor operator*(double x, const Tensor& t) {
    
  Tensor res(t.get_mp(), t.get_valence(), t.get_index_type(), 
		    t.get_triad()) ;

  for (int i=0 ; i<res.get_n_comp() ; i++) {
    Itbl ind (res.indices(i)) ;
    res.set(ind) = x*t(ind) ;
  }

  return res ; 

}


Tensor operator* (const Tensor& t, double x) {
    return x * t ;
}

Tensor operator*(int m, const Tensor& t) {
    return double(m) * t ; 
}


Tensor operator* (const Tensor& t, int m) {
    return double(m) * t ;
}


			//**********//
			// DIVISION //
			//**********//

Tensor operator/(const Tensor& t1, const Scalar& s2) {
    
    // Protections
    assert(s2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mp() == s2.get_mp()) ;

    // Cas particuliers
    if (s2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Tensor / Scalar !" << endl ;
	abort() ; 
    }

    if (s2.get_etat() == ETATUN) return t1 ;

    Tensor res(t1.get_mp(), t1.get_valence(), t1.get_index_type(), 
		t1.get_triad()) ;

    for (int i=0 ; i<res.get_n_comp() ; i++) {
	Itbl ind (res.indices(i)) ;
	res.set(ind) = t1(ind) / s2 ;	    // Scalar / Scalar
    }
    return res ;

}


Tensor operator/ (const Tensor& t, double x) {

    if ( x == double(0) ) {
	cout << "Division by 0 in Tensor / double !" << endl ;
	abort() ;
    }

    if (x == double(1)) 
      return t ;
    else {
	Tensor res(t.get_mp(), t.get_valence(), t.get_index_type(), 
		    t.get_triad()) ;

	for (int i=0 ; i<res.get_n_comp() ; i++) {
	  Itbl ind (res.indices(i)) ;
	  res.set(ind) = t(ind) / x ;	    // Scalar / double
	}
	return res ; 
    }

}

Tensor operator/ (const Tensor& t, int m) {

    return t / double(m) ; 
}







}
