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

char tensor_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/26 14:33:53  j_novak
 * Arithmetic functions for the class Tensor
 *
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
#include "tensor.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

Tensor operator+(const Tensor & t) {

    return t ; 

}

Tensor operator-(const Tensor & t) {
    
  Tensor res(*(t.get_mp()), t.get_valence(), t.get_index_type(), 
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

    Tensor res(*(t1.get_mp()), t1.get_valence(), t1.get_index_type(), 
			t1.get_triad()) ;

    for (int i=0 ; i<res.get_n_comp() ; i++) {
      Itbl ind (res.indices(i)) ;
      res.set(ind) = t1(ind) + t2(ind) ;
    }
    return res ;

}

			//**************//
			// SOUSTRACTION //
			//**************//

Tensor operator-(const Tensor & t1, const Tensor & t2) {

    return (t1 + (-t2)) ;

}

			//****************//
			// MULTIPLICATION //
			//****************//



Tensor operator*(double x, const Tensor& t) {
    
  Tensor res(*(t.get_mp()), t.get_valence(), t.get_index_type(), 
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

Tensor operator/ (const Tensor& t1, const Scalar& s2) {
    
    // Protections
    assert(s2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mp() == s2.get_mp()) ;

    // Cas particuliers
    if (s2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Tensor / Tensor !" << endl ;
	abort() ; 
    }
    // Cas general
    assert(s2.get_etat() == ETATQCQ) ;  // sinon...

    Tensor res(*(t1.get_mp()), t1.get_valence(), t1.get_index_type(), 
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
	Tensor res(*(t.get_mp()), t.get_valence(), t.get_index_type(), 
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







