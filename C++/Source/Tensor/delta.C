/*
 *  Methods of class Delta
 *
 *   (see file tensor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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


char delta_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/29 13:48:18  j_novak
 * New class Delta.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <assert.h>

// Headers Lorene
#include "tensor.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Delta::Delta(const Map& map, const Itbl& tipe,const Base_vect& triad_i) 
		: Tensor(map, 3, tipe, 18, triad_i) {

}

// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Delta::Delta(const Map& map, int tipe, const Base_vect& triad_i)  
  : Tensor(map, 3, tipe, 18, triad_i){

}

// Copy constructor
// ----------------
Delta::Delta (const Delta& source) : 
  Tensor (*source.mp, 3, source.type_indice, 18, *(source.triad)) {
    
    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.position(indices(i)) ;
	cmp[i] = new Scalar (*source.cmp[place_source]) ;
    }
}   


// Constructor from a Tensor
// --------------------------
Delta::Delta (const Tensor& source) :
  Tensor (*source.mp, 3, source.type_indice, 18, *(source.triad)) {
	
    assert (source.valence == 3) ;

    for (int i=0 ; i<n_comp ; i++) {
	int place_source = source.position(indices(i)) ;
	cmp[i] = new Scalar (*source.cmp[place_source]) ;
    }
	    
}   

	
// Constructor from a file
// -----------------------
Delta::Delta(const Map& map, const Base_vect& triad_i, FILE* fd)
			: Tensor(map, triad_i, fd) {
	
	assert (valence == 3) ;
	assert (n_comp == 18) ;
}

			//--------------//
			//  Destructor  //
			//--------------//

Delta::~Delta() {}




	
int Delta::position (const Itbl& idx) const {
    
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == 3) ;
    for (int i=0 ; i<3 ; i++)
	assert ((idx(i) >= 1) && (idx(i) <= 3)) ;
	
     // Gestion des deux derniers indices :
    int last = idx(2) ;
    int first = idx(1) ;
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
    
    return 6*(idx(0)-1) + place_fin ;
}

Itbl Delta::indices (int place) const {
    Itbl res(3) ;
    res.set_etat_qcq() ;
    assert ((place>=0) && (place<18)) ;
    
    res.set(0) = place / 6 + 1 ;
    if (place<3) {
	res.set(1) = 1 ;
	res.set(2) = place+1 ;
	}
    
    if ((place>2) && (place<5)) {
	res.set(1) = 2 ;
	res.set(2) = place - 1 ;
	}
    
    if (place == 5) {
	res.set(1) = 3 ;
	res.set(2) = 3 ;
	}
 
    return res ;
}
	
void Delta::operator= (const Tensor& t) {
    
    assert (t.get_valence() == 3) ;
    
    triad = t.triad ; 
    
    for (int i=0 ; i<3 ; i++)
      assert (type_indice(i) == t.type_indice(i)) ;
    
    for (int i=0 ; i<18 ; i++) {
      int place_t = t.position(indices(i)) ;
      *cmp[i] = *t.cmp[place_t] ;
    }
}

