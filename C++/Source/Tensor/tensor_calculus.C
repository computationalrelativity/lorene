/*
 *  Methods of class Tensor for tensor calculus
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


char tensor_calculus_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
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


Tensor Tensor::contract(int ind_1, int ind_2) const {
    
    // Les verifications :
    assert ((ind_1 >= 0) && (ind_1 < valence)) ;
    assert ((ind_2 >= 0) && (ind_2 < valence)) ;
    assert (type_indice(ind_1) != type_indice(ind_2))  ;
 
    // On veut ind_1 < ind_2 :
    if (ind_1 > ind_2) {
		int auxi = ind_2 ;
		ind_2 = ind_1 ;
		ind_1 = auxi ;
    }
    
    // On construit le resultat :
    int val_res = valence - 2 ;
   
    Itbl tipe(val_res) ;
    tipe.set_etat_qcq() ;
	
    for (int i=0 ; i<ind_1 ; i++)
		tipe.set(i) = type_indice(i) ;
    for (int i=ind_1 ; i<ind_2-1 ; i++)
		tipe.set(i) = type_indice(i+1) ;
    for (int i = ind_2-1 ; i<val_res ; i++)
		tipe.set(i) = type_indice(i+2) ;
	
    Tensor res(*mp, val_res, tipe, triad) ; 
	
    Scalar work(*mp) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_source(valence) ;
    jeux_indice_source.set_etat_qcq() ;
	
    for (int i=0 ; i<res.get_n_comp() ; i++) {

		Itbl jeux_indice_res(res.indices(i)) ;
		
		for (int j=0 ; j<ind_1 ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j) ;
		for (int j=ind_1+1 ; j<ind_2 ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j-1) ;
		for (int j=ind_2+1 ; j<valence ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j-2) ;
	    
		work.set_etat_zero() ;
		for (int j=1 ; j<=3 ; j++) {
	    	jeux_indice_source.set(ind_1) = j ;
	    	jeux_indice_source.set(ind_2) = j ;
	    	work = work + (*cmp[position(jeux_indice_source)]) ;
	    }
	    
		res.set(jeux_indice_res) = work ;
	}
	
    return res ;
}



				//----------------------//
				//  Index manipulation	//
				//----------------------//


Tensor Tensor::up(int , const Metric& ) const {
	
		cout << " Tensor::up :  not ready yet ! " << endl ; 
		abort() ; 
		return *this ; 
} 

















 
