/*
 *  Function external to class Tensor for tensor calculus
 *
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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


char tensor_calculus_ext_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/11 16:47:10  e_gourgoulhon
 * Suppressed the call to Ibtl::set_etat_qcq() after the construction
 * of the Itbl's, thanks to the new property of the Itbl class.
 *
 * Revision 1.1  2003/10/06 15:13:38  e_gourgoulhon
 * Tensor contraction.
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

				//------------------//
				//   Contraction	//
				//------------------//


Tensor contract(const Tensor& t1, int ind1, const Tensor& t2, int ind2) {
    
	int val1 = t1.get_valence() ; 
	int val2 = t2.get_valence() ; 

    // Verifs :
    assert((ind1>=0) && (ind1<val1)) ;
    assert((ind2>=0) && (ind2<val2)) ;
    assert(t1.get_mp() == t2.get_mp()) ;
    
    // Contraction possible ?
    if ( (val1 != 0) && (val2 != 0) ) {
	    assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    assert (t1.get_index_type(ind1) != t2.get_index_type(ind2)) ;
    
    int val_res = val1 + val2 - 2;
	
    Itbl tipe(val_res) ;

    for (int i=0 ; i<ind1 ; i++)
		tipe.set(i) = t1.get_index_type(i) ;
    for (int i=ind1 ; i<val1-1 ; i++)
		tipe.set(i) = t1.get_index_type(i+1) ;
    for (int i=val1-1 ; i<val1+ind2-1 ; i++)
		tipe.set(i) = t2.get_index_type(i-val1+1) ;
    for (int i = val1+ind2-1 ; i<val_res ; i++)
		tipe.set(i) = t2.get_index_type(i-val1+2) ;
	
    const Base_vect* triad_res = (val_res == 0) ? 0x0 : t1.get_triad() ; 

    Tensor res(t1.get_mp(), val_res, tipe, triad_res) ;
	
    Scalar work(t1.get_mp()) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_t1(val1) ;
    Itbl jeux_indice_t2(val2) ;
    
    for (int i=0 ; i<res.get_n_comp() ; i++) {
	
		Itbl jeux_indice_res(res.indices(i)) ;
		
		for (int j=0 ; j<ind1 ; j++)
	    	jeux_indice_t1.set(j) = jeux_indice_res(j) ;
			
		for (int j=ind1+1 ; j<val1 ; j++)
	    	jeux_indice_t1.set(j) = jeux_indice_res(j-1) ;

		for (int j=0 ; j<ind2 ; j++)
	    	jeux_indice_t2.set(j) = jeux_indice_res(val1+j-1) ;

		for (int j=ind2+1 ; j<val2 ; j++)
	    	jeux_indice_t2.set(j) = jeux_indice_res(val1+j-2) ;
	
		work.set_etat_zero() ;
		for (int j=1 ; j<=3 ; j++) {
	    	jeux_indice_t1.set(ind1) = j ;
	    	jeux_indice_t2.set(ind2) = j ;
	    	work = work + t1(jeux_indice_t1) * t2(jeux_indice_t2) ;
	    }
	    
		res.set(jeux_indice_res) = work ;
	}
	
    return res ;
}







