/*
 *  Tensor calculus for class Tensor_sym
 *
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


char tensor_sym_calculus_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/01/30 12:44:53  e_gourgoulhon
 * Added Tensor_sym operator*(const Tensor_sym&, const Tensor_sym& ).
 *
 * Revision 1.1  2004/01/08 09:22:40  e_gourgoulhon
 * First version.
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

// Tensorial product
//------------------

Tensor_sym operator*(const Tensor_sym& t1, const Tensor& t2) {
   
    assert (t1.mp == t2.mp) ;
    
    int val_res = t1.valence + t2.valence ;
     
    Itbl tipe(val_res) ;
  
    for (int i=0 ; i<t1.valence ; i++)
	    tipe.set(i) = t1.type_indice(i) ;
    for (int i=0 ; i<t2.valence ; i++)
	    tipe.set(i+t1.valence) = t2.type_indice(i) ;
    
    
    const Base_vect* triad_res = t1.get_triad() ; 
    
    if ( t2.valence != 0 ) {
	    assert ( *(t2.get_triad()) == *triad_res ) ;
    }
    
    Tensor_sym res(*t1.mp, val_res, tipe, *triad_res, t1.sym_index1(),
                    t1.sym_index2()) ;
    
    Itbl jeux_indice_t1(t1.valence) ;
    Itbl jeux_indice_t2(t2.valence) ;
        
    for (int i=0 ; i<res.n_comp ; i++) {
	    Itbl jeux_indice_res(res.indices(i)) ;
	    for (int j=0 ; j<t1.valence ; j++)
	        jeux_indice_t1.set(j) = jeux_indice_res(j) ;
	    for (int j=0 ; j<t2.valence ; j++)
	        jeux_indice_t2.set(j) = jeux_indice_res(j+t1.valence) ;
	
	    res.set(jeux_indice_res) = t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
    }
    
    return res ;
}


Tensor_sym operator*(const Tensor& t1, const Tensor_sym& t2) {
   
    assert (t1.mp == t2.mp) ;
    
    int val_res = t1.valence + t2.valence ;
     
    Itbl tipe(val_res) ;
  
    for (int i=0 ; i<t1.valence ; i++)
	    tipe.set(i) = t1.type_indice(i) ;
    for (int i=0 ; i<t2.valence ; i++)
	    tipe.set(i+t1.valence) = t2.type_indice(i) ;
    
    
    const Base_vect* triad_res = t2.get_triad() ; 
    
    if ( t1.valence != 0 ) {
	    assert ( *(t1.get_triad()) == *triad_res ) ;
    }
    
    int ids1 = t2.sym_index1() + t1.valence ; // symmetry index 1 of the result
    int ids2 = t2.sym_index2() + t1.valence ; // symmetry index 2 of the result

    Tensor_sym res(*t2.mp, val_res, tipe, *triad_res, ids1, ids2) ;
    
    Itbl jeux_indice_t1(t1.valence) ;
    Itbl jeux_indice_t2(t2.valence) ;
        
    for (int i=0 ; i<res.n_comp ; i++) {
	    Itbl jeux_indice_res(res.indices(i)) ;
	    for (int j=0 ; j<t1.valence ; j++)
	        jeux_indice_t1.set(j) = jeux_indice_res(j) ;
	    for (int j=0 ; j<t2.valence ; j++)
	        jeux_indice_t2.set(j) = jeux_indice_res(j+t1.valence) ;
	
	    res.set(jeux_indice_res) = t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
    }
    
    return res ;
}



Tensor_sym operator*(const Tensor_sym& t1, const Tensor_sym& t2) {
   
    assert (t1.mp == t2.mp) ;
    
    int val_res = t1.valence + t2.valence ;
     
    Itbl tipe(val_res) ;
  
    for (int i=0 ; i<t1.valence ; i++)
	    tipe.set(i) = t1.type_indice(i) ;
    for (int i=0 ; i<t2.valence ; i++)
	    tipe.set(i+t1.valence) = t2.type_indice(i) ;
    
    
    const Base_vect* triad_res = t1.get_triad() ; 
    
	assert ( *(t2.get_triad()) == *triad_res ) ;
    
    Tensor_sym res(*t1.mp, val_res, tipe, *triad_res, t1.sym_index1(),
                    t1.sym_index2()) ;
    
    Itbl jeux_indice_t1(t1.valence) ;
    Itbl jeux_indice_t2(t2.valence) ;
        
    for (int i=0 ; i<res.n_comp ; i++) {
	    Itbl jeux_indice_res(res.indices(i)) ;
	    for (int j=0 ; j<t1.valence ; j++)
	        jeux_indice_t1.set(j) = jeux_indice_res(j) ;
	    for (int j=0 ; j<t2.valence ; j++)
	        jeux_indice_t2.set(j) = jeux_indice_res(j+t1.valence) ;
	
	    res.set(jeux_indice_res) = t1(jeux_indice_t1)*t2(jeux_indice_t2) ;
    }
    
    return res ;
}














 
