/*
 *  Arithmetics Tensor_sym
 *
 *   (see file tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (for preceding class Tenseur)
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


char tensor_sym_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:20  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2004/01/08 09:22:40  e_gourgoulhon
 * First version.
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
Tensor_sym operator+(const Tensor_sym & t) {

    return t ; 

}


Tensor_sym operator-(const Tensor_sym & tt) {
    
  Tensor_sym res(tt.get_mp(), tt.get_valence(), tt.get_index_type(), 
		    *(tt.get_triad()), tt.sym_index1(), tt.sym_index2()) ;

  for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
    Itbl ind = res.indices(ic) ;    
    res.set(ind) = -tt(ind) ;
  }
  return res ;

}

			//**********//
			// ADDITION //
			//**********//


Tensor_sym operator+(const Tensor_sym& t1, const Tensor_sym& t2) {
    
    assert (t1.get_valence() == t2.get_valence()) ;
    assert (t1.get_mp() == t2.get_mp()) ;
    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    
    for (int id=0 ; id<t1.get_valence() ; id++)
	assert(t1.get_index_type(id) == t2.get_index_type(id)) ;
        
    int ids1 = t1.sym_index1() ; 
    int ids2 = t1.sym_index2() ;
    
    assert(t2.sym_index1() == ids1) ; 
    assert(t2.sym_index2() == ids2) ; 

    Tensor_sym res(t1.get_mp(), t1.get_valence(), t1.get_index_type(), 
		   *(t1.get_triad()), ids1, ids2) ;

    for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
      Itbl ind = res.indices(ic) ;
      res.set(ind) = t1(ind) + t2(ind) ;
    }
    return res ;

}

			//**************//
			// SOUSTRACTION //
			//**************//


Tensor_sym operator-(const Tensor_sym& t1, const Tensor_sym& t2) {
    
    assert (t1.get_valence() == t2.get_valence()) ;
    assert (t1.get_mp() == t2.get_mp()) ;
    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    
    for (int id=0 ; id<t1.get_valence() ; id++)
	assert(t1.get_index_type(id) == t2.get_index_type(id)) ;
        
    int ids1 = t1.sym_index1() ; 
    int ids2 = t1.sym_index2() ;
    
    assert(t2.sym_index1() == ids1) ; 
    assert(t2.sym_index2() == ids2) ; 

    Tensor_sym res(t1.get_mp(), t1.get_valence(), t1.get_index_type(), 
		   *(t1.get_triad()), ids1, ids2) ;

    for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
      Itbl ind = res.indices(ic) ;
      res.set(ind) = t1(ind) - t2(ind) ;
    }
    return res ;

}



			//****************//
			// MULTIPLICATION //
			//****************//


Tensor_sym operator*(const Scalar& t1, const Tensor_sym& t2) {
   
    assert (&(t1.get_mp()) == &(t2.get_mp())) ;
        
    if (t1.get_etat() == ETATUN) return t2 ;
    
    Tensor_sym res(t2.get_mp(), t2.get_valence(), t2.get_index_type(), 
               *(t2.get_triad()), t2.sym_index1(), t2.sym_index2()) ;
    
    for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
        Itbl ind = res.indices(ic) ;
        res.set(ind) = t1 * t2(ind) ;
    }
    
    return res ;
}


Tensor_sym operator*(const Tensor_sym& t2, const Scalar& t1) {
       
    return t1*t2 ;
}



Tensor_sym operator*(double x, const Tensor_sym& tt) {
    
    Tensor_sym res(tt.get_mp(), tt.get_valence(), tt.get_index_type(), 
		    *(tt.get_triad()), tt.sym_index1(), tt.sym_index2()) ;

    for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
        Itbl ind = res.indices(ic) ;
        res.set(ind) = x * tt(ind) ;
    }

    return res ; 

}


Tensor_sym operator*(const Tensor_sym& t, double x) {
    return x * t ;
}


Tensor_sym operator*(int m, const Tensor_sym& t) {
    return double(m) * t ; 
}


Tensor_sym operator*(const Tensor_sym& t, int m) {
    return double(m) * t ;
}


			//**********//
			// DIVISION //
			//**********//

Tensor_sym operator/(const Tensor_sym& t1, const Scalar& s2) {
    
    // Protections
    assert(s2.get_etat() != ETATNONDEF) ;
    assert(t1.get_mp() == s2.get_mp()) ;

    // Cas particuliers
    if (s2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Tensor_sym / Scalar !" << endl ;
	abort() ; 
    }

    if (s2.get_etat() == ETATUN) return t1 ;

    Tensor_sym res(t1.get_mp(), t1.get_valence(), t1.get_index_type(), 
		    *(t1.get_triad()), t1.sym_index1(), t1.sym_index2()) ;

    for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
        Itbl ind = res.indices(ic) ;
        res.set(ind) = t1(ind) / s2 ;
    }
    
    return res ;

}


Tensor_sym operator/(const Tensor_sym& tt, double x) {

    if ( x == double(0) ) {
	cout << "Division by 0 in Tensor_sym / double !" << endl ;
	abort() ;
    }

    if ( x == double(1) ) return tt ;
    else {
        Tensor_sym res(tt.get_mp(), tt.get_valence(), tt.get_index_type(), 
		    *(tt.get_triad()), tt.sym_index1(), tt.sym_index2()) ;

        for (int ic=0 ; ic<res.get_n_comp() ; ic++) {
            Itbl ind = res.indices(ic) ;
            res.set(ind) = tt(ind) / x ;
        }

	return res ; 
    }

}


Tensor_sym operator/(const Tensor_sym& t, int m) {

    return t / double(m) ; 
}


}
