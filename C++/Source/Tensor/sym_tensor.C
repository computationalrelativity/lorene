/*
 *  Methods of class Sym_tensor
 *
 *   (see file tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (Cmp version)
 *   Copyright (c) 2000-2001 Eric Gourgoulhon (Cmp version)
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


char sym_tensor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2003/10/11 16:47:10  e_gourgoulhon
 * Suppressed the call to Ibtl::set_etat_qcq() after the construction
 * of the Itbl's, thanks to the new property of the Itbl class.
 *
 * Revision 1.7  2003/10/07 09:56:59  j_novak
 * method Sym_tensor::inverse() implemented (but not tested!)
 *
 * Revision 1.6  2003/10/06 13:58:48  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.5  2003/10/03 11:21:48  j_novak
 * More methods for the class Metric
 *
 * Revision 1.4  2003/10/02 15:45:51  j_novak
 * New class Metric
 *
 * Revision 1.3  2003/10/01 15:39:43  e_gourgoulhon
 * Added assert to insure that both indices have the same type.
 *
 * Revision 1.2  2003/09/26 08:05:31  j_novak
 * New class Vector.
 *
 * Revision 1.1  2003/09/25 13:37:40  j_novak
 * Symmetric tensors of valence 2 are now implemented (not tested yet).
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

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Sym_tensor::Sym_tensor(const Map& map, const Itbl& tipe,const Base_vect& triad_i) 
		: Tensor(map, 2, tipe, 6, triad_i) {
		
		assert(tipe(0) == tipe(1)) ; 
		set_der_0x0() ;

}

// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Sym_tensor::Sym_tensor(const Map& map, int tipe, const Base_vect& triad_i)  
  : Tensor(map, 2, tipe, 6, triad_i){

		set_der_0x0() ;
}

// Copy constructor
// ----------------
Sym_tensor::Sym_tensor (const Sym_tensor& source) : 
  Tensor (*source.mp, 2, source.type_indice, 6, *(source.triad)) {
    
    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.position(indices(i)) ;
	cmp[i] = new Scalar (*source.cmp[place_source]) ;
    }
		set_der_0x0() ;
}   


// Constructor from a Tensor
// --------------------------
Sym_tensor::Sym_tensor (const Tensor& source) :
  Tensor (*source.mp, 2, source.type_indice, 6, *(source.triad)) {
	
    assert (source.valence == 2) ;
	assert(source.type_indice(0) == source.type_indice(1)) ; 
	

    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.position(indices(i)) ;
	cmp[i] = new Scalar (*source.cmp[place_source]) ;
    }
	    
		set_der_0x0() ;
}   

	
// Constructor from a file
// -----------------------
Sym_tensor::Sym_tensor(const Map& map, const Base_vect& triad_i, FILE* fd)
			: Tensor(map, triad_i, fd) {
	
	assert (valence == 2) ;
	assert (n_comp == 6) ;
	set_der_0x0() ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Sym_tensor::~Sym_tensor() {

  del_deriv() ;

}



			//-------------------//
			// Memory managment  //
			//-------------------//

void Sym_tensor::del_deriv() const {

  set_der_0x0() ;
  Tensor::del_deriv() ;

}

void Sym_tensor::set_der_0x0() const {

}


	
int Sym_tensor::position (const Itbl& idx) const {
    
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == 2) ;
    for (int i=0 ; i<2 ; i++)
	assert ((idx(i) >= 1) && (idx(i) <= 3)) ;
	
     // Gestion des deux indices :
    int last = idx(1) ;
    int first = idx(0) ;
    if (last < first) {
	int auxi = last ;
	last = first ;
	first = auxi ;
    }
    
    int place_fin ;
    switch (first) {
    case 1 : {
      place_fin = last - 1;
      break ;
    }
    case 2 : {
      place_fin = 1+last ;
      break ;
    }
    case 3 : {
      place_fin = 5 ;
      break ;
    }
    default : {
      abort() ;
    }
    }
    
    return place_fin ;
}

Itbl Sym_tensor::indices (int place) const {
    Itbl res(2) ;
    assert ((place>=0) && (place<6)) ;
    
    if (place<3) {
		res.set(0) = 1 ;
		res.set(1) = place+1 ;
	}
    
    if ((place>2) && (place<5)) {
		res.set(0) = 2 ;
		res.set(1) = place - 1 ;
	}
    
    if (place == 5) {
		res.set(0) = 3 ;
		res.set(1) = 3 ;
	}
 
    return res ;
}
	
void Sym_tensor::operator= (const Tensor& t) {
    
    assert (t.get_valence() == 2) ;
    
    triad = t.triad ; 
    
    for (int i=0 ; i<2 ; i++)
      assert (type_indice(i) == t.type_indice(i)) ;
    
    for (int i=0 ; i<6 ; i++) {
      int place_t = t.position(indices(i)) ;
      *cmp[i] = *t.cmp[place_t] ;
    }

    del_deriv() ;
}

Sym_tensor* Sym_tensor::inverse() const {

  //Le resultat :
  Sym_tensor* res = 
    new Sym_tensor(*mp, -type_indice(0), *triad) ;
    
  // le determinant :
  Scalar determ1(*mp) ;
  determ1 = double(1)/
    (operator()(1, 1)*operator()(2, 2)*operator()(3, 3) 
     + operator()(1, 2)*operator()(2, 3)*operator()(1, 3)
     + operator()(1, 3)*operator()(1, 2)*operator()(2, 3) 
     - operator()(1, 3)*operator()(2, 2)*operator()(1, 3)
     - operator()(2, 3)*operator()(2, 3)*operator()(1, 1) 
     - operator()(3, 3)*operator()(1, 2)*operator()(1, 2) ) ;
    
  int sgn ;	// Le signe du co-facteur ...
  int l_up, l_down, c_left, c_right ;	    // Coordonnees du cofacteur :
    
  Scalar cofacteur(*mp) ;
    
  for (int i=1 ; i<=3 ; i++) {
    sgn = 1 ;
    for (int j=i ; j<=3 ; j++) {
	    
      switch (j) {
		
      case 1 : {
	c_left = 2 ;
	c_right = 3 ;
	break ;
      }
      case 2 : {
	c_left = 1 ;
	c_right = 3 ;
	break ;
      }
      default : {
	c_left = 1 ;
	c_right = 2 ;
	break ;
      }
      }
	    
      switch (i) {
		
      case 1 : {
	l_up = 2 ;
	l_down = 3 ;
	break ;
      }
      case 2 : {
	l_up = 1 ;
	l_down = 3 ;
	break ;
      }
      default : {
	l_up = 1 ;
	l_down = 2 ;
	break ;
      } 
      }
	    
      cofacteur = sgn*(operator()(l_up, c_left)*operator()(l_down, c_right)-
		       operator()(l_up, c_right)*operator()(l_down, c_left))*determ1 ;
	    
      res->set(i, j) = cofacteur ;
      sgn *= -1 ;
    }
  }
  return res ;

}

