/*
 *  Methods of class Tensor_sym
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


char tensor_sym_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/01/04 20:51:45  e_gourgoulhon
 * New class to deal with general tensors which are symmetric with
 * respect to two of their indices.
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
#include "utilitaires.h"


			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
Tensor_sym::Tensor_sym(const Map& map, int val, const Itbl& tipe, 
		 const Base_vect& triad_i, int index_sym1, int index_sym2) 
            : Tensor(map, val, tipe, 6*int(pow(3.,val-2)), triad_i),
              id_sym1(index_sym1),    
              id_sym2(index_sym2) {    
		
    // Des verifs :
    assert( valence >= 2 ) ;
    assert( (id_sym1 >=0) && (id_sym1 < valence) ) ;  
    assert( (id_sym2 >=0) && (id_sym2 < valence) ) ;  
    assert( id_sym1 != id_sym2 ) ; 

    // The symmetry indices must be of same type: 
    assert( tipe(id_sym1) == tipe(id_sym2) ) ; 
    
    // Possible re-ordering of the symmetry indices
    if (id_sym1 > id_sym2) {
        int tmp = id_sym1 ; 
        id_sym1 = id_sym2 ; 
        id_sym2 = tmp ; 
    }
    	
}



// Standard constructor when all the indices are of the same type
// --------------------------------------------------------------
Tensor_sym::Tensor_sym(const Map& map, int val, int tipe, 
		 const Base_vect& triad_i, int index_sym1, int index_sym2) 
            : Tensor(map, val, tipe, 6*int(pow(3.,val-2)), triad_i),
              id_sym1(index_sym1),    
              id_sym2(index_sym2) {    
		
    // Des verifs :
    assert( valence >= 2 ) ;
    assert( (id_sym1 >=0) && (id_sym1 < valence) ) ;  
    assert( (id_sym2 >=0) && (id_sym2 < valence) ) ;  
    assert( id_sym1 != id_sym2 ) ; 

    // Possible re-ordering of the symmetry indices
    if (id_sym1 > id_sym2) {
        int tmp = id_sym1 ; 
        id_sym1 = id_sym2 ; 
        id_sym2 = tmp ; 
    }
    	
}

// Constructor for a valence 3 symmetric tensor
// --------------------------------------------
Tensor_sym::Tensor_sym(const Map& map, int tipe0, int tipe1, int tipe2, 
                   const Base_vect& triad_i,
                   int index_sym1, int index_sym2) 
            : Tensor(map, 3, tipe0, 18, triad_i),
              id_sym1(index_sym1),    
              id_sym2(index_sym2) {    
    
    assert( (tipe0==COV) || (tipe0==CON) ) ;		
    assert( (tipe1==COV) || (tipe1==CON) ) ;		
    assert( (tipe2==COV) || (tipe2==CON) ) ;		

    type_indice.set(1) = tipe1 ; 
    type_indice.set(2) = tipe2 ; 

    assert( (id_sym1 >=0) && (id_sym1 < 3) ) ;  
    assert( (id_sym2 >=0) && (id_sym2 < 3) ) ;  
    assert( id_sym1 != id_sym2 ) ; 
    assert( type_indice(id_sym1) == type_indice(id_sym2) ) ; 

    // Possible re-ordering of the symmetry indices
    if (id_sym1 > id_sym2) {
        int tmp = id_sym1 ; 
        id_sym1 = id_sym2 ; 
        id_sym2 = tmp ; 
    }
    	
}



// Copy constructor
// ----------------
Tensor_sym::Tensor_sym(const Tensor_sym& source) 
            : Tensor(*source.mp, source.valence, source.type_indice, 
                     6*int(pow(3.,source.valence-2)) , *(source.triad)),
              id_sym1(source.id_sym1),    
              id_sym2(source.id_sym2) {    
                     
    for (int i=0 ; i<n_comp ; i++) {

        int posi = source.position(indices(i)) ;  // in case source belongs to
                                                  // a derived class of 
                                                  // Tensor_sym with a different
                                                  // storage of components 
	*(cmp[i]) = *(source.cmp[posi]) ;
    }
}   



	

// Constructor from a file
// -----------------------
Tensor_sym::Tensor_sym(const Map& map, const Base_vect& triad_i, FILE* fd)
			: Tensor(map, triad_i, fd) {
	
    fread_be(&id_sym1, sizeof(int), 1, fd) ;
    fread_be(&id_sym2, sizeof(int), 1, fd) ;
    
    assert( type_indice(id_sym1) == type_indice(id_sym2) ) ; 
    
}


			//--------------//
			//  Destructor  //
			//--------------//

Tensor_sym::~Tensor_sym() {

}

			//--------------//
			//  Assignment  //
			//--------------//


void Tensor_sym::operator=(const Tensor_sym& tt) {
    
    assert (valence == tt.valence) ;

    triad = tt.triad ; 
    id_sym1 = tt.id_sym1 ; 
    id_sym2 = tt.id_sym2 ; 

    for (int id=0 ; id<valence ; id++)
      assert(tt.type_indice(id) == type_indice(id)) ;
	
    for (int ic=0 ; ic<n_comp ; ic++) {
        int posi = tt.position(indices(ic)) ;
        *cmp[ic] = *(tt.cmp[posi]) ;
    }

    del_deriv() ;

}

void Tensor_sym::operator=(const Tensor& tt) {
    
    assert (valence == tt.get_valence()) ;

    triad = tt.get_triad() ; 

    for (int id=0 ; id<valence ; id++)
      assert(tt.type_indice(id) == type_indice(id)) ;

    // The symmetry indices must be of same type: 
    assert( tt.type_indice(id_sym1) == tt.type_indice(id_sym2) ) ; 

	
    for (int ic=0 ; ic<n_comp ; ic++) {
        int posi = tt.position(indices(ic)) ;
        *cmp[ic] = *(tt.cmp[posi]) ;
    }

    del_deriv() ;

}


			//--------------//
			//   Accessor   //
			//--------------//

int Tensor_sym::position(const Itbl& idx) const {
    
    // Protections:
    assert (idx.get_ndim() == 1) ;
    assert (idx.get_dim(0) == valence) ;
    for (int i=0 ; i<valence ; i++) {
	assert( (idx(i)>=1) && (idx(i)<=3) ) ;
    }
    
    // The two symmetric indices are moved to the end --> new index array idx0
    Itbl idx0(valence) ; 
    if (valence > 2) {
        for (int id=0 ; id<id_sym1; id++) {
            idx0.set(id) = idx(id) ; 
        }
        for (int id=id_sym1; id<id_sym2-1; id++) {
            idx0.set(id) = idx(id+1) ;
        }
        for (int id=id_sym2-1; id<valence-2; id++) {
            idx0.set(id) = idx(id+2) ; 
        }
        idx0.set(valence-2) = idx(id_sym1) ;  //## not used
        idx0.set(valence-1) = idx(id_sym2) ;  //## in what follows
    }
    
    // Values of the symmetric indices:
    int is1 = idx(id_sym1) ; 
    int is2 = idx(id_sym2) ; 
    
    // Reordering to ensure is1 <= is2 :
    if (is2 < is1) {
        int aux = is1 ; 
        is1 = is2 ; 
        is2 = aux ; 
    }

    // Position in the cmp array :
    int pos = 0 ; 
    for (int id=0 ; id<valence-2 ; id++) {
        pos = 3 * pos + idx0(id) - 1 ;  // all the values of each non symmetric
                                        // index occupy 3 "boxes"
    }
    
    pos = 6 * pos  ;   // all the values of the two symmetric
                       // indices occupy 6 "boxes"
    switch (is1) {
        case 1 : {
            pos += is2 - 1 ;     // (1,1), (1,2) and (1,3) stored respectively
            break ;              // in relative position 0, 1 and 2
        }        
        case 2 : {
            pos += is2 + 1 ;     // (2,2) and (2,3) stored respectively
            break ;              // in relative position 3 and  4
        }
        case 3 : {
            pos += 5 ;           // (3,3) stored in relative position 5
            break ; 
        }
    }
    
    return pos ;
}



Itbl Tensor_sym::indices(int place) const {

    assert( (place>=0) && (place<n_comp) ) ;
    
    // Index set with the two symmetric indices at the end:

    Itbl idx0(valence) ; 
    
    int reste = div(place, 6).rem ;
    place = int((place-reste)/6) ;
    
    if (reste<3) {
        idx0.set(valence-2) = 1 ;
	idx0.set(valence-1) = reste + 1 ;
    }
    
    if ( (reste>2) && (reste<5) ) {
	idx0.set(valence-2) = 2 ;
	idx0.set(valence-1) = reste - 1 ;
    }
    
    if (reste == 5) {
        idx0.set(valence-2) = 3 ;
        idx0.set(valence-1) = 3 ;
    }
    
    // The output is ready in the case of a valence 2 tensor: 
    if (valence == 2) return idx0 ; 
    
    for (int id=valence-3 ; id>=0 ; id--) {
	int ind = div(place, 3).rem ;
	place = int((place-ind)/3) ;
	idx0.set(id) = ind + 1 ; 
    }
	
    // Reorganization of the index set to put the two symmetric indices at
    // their correct positions:
    
    Itbl idx(valence) ; 
    
    for (int id=0 ; id<id_sym1; id++) {
        idx.set(id) = idx0(id) ; 
    }
    idx.set(id_sym1) = idx0(valence-2) ; 
    
    for (int id=id_sym1+1; id<id_sym2; id++) {
        idx.set(id) = idx0(id-1) ;
    }
    idx.set(id_sym2) = idx0(valence-1) ; 
    
    for (int id=id_sym2+1; id<valence; id++) {
        idx.set(id) = idx0(id-2) ; 
    }

    return idx ;
}
	

			//--------------//
			//    Outputs   //
			//--------------//

void Tensor_sym::sauve(FILE* fd) const {

    Tensor::sauve(fd) ; 
    
    fwrite_be(&id_sym1, sizeof(int), 1, fd) ;   
    fwrite_be(&id_sym2, sizeof(int), 1, fd) ;   

}









