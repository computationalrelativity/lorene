/*
 *  Function external to class Tensor for tensor calculus
 *
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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
 * Revision 1.6  2004/01/23 08:00:16  e_gourgoulhon
 * Minor modifs. in output of methods min, max, maxabs, diffrel to
 * better handle the display in the scalar case.
 *
 * Revision 1.5  2004/01/15 10:59:53  f_limousin
 * Added method contract_desal for the contraction of two tensors with desaliasing
 *
 * Revision 1.4  2004/01/14 11:38:32  f_limousin
 * Added method contract for one tensor
 *
 * Revision 1.3  2003/11/05 15:29:36  e_gourgoulhon
 *  Added declaration of externa functions max, min, maxabs,
 * diffrel and diffrelmax.
 *
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
				//   Contraction    //
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


Tensor contract_desal(const Tensor& t1, int ind1, const Tensor& t2, int ind2) {
    
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
	    	work = work + t1(jeux_indice_t1) % t2(jeux_indice_t2) ;
	    }
	    
		res.set(jeux_indice_res) = work ;
	}
	
    return res ;
}


Tensor contract(const Tensor& source, int ind_1, int ind_2) {
    
    int val = source.get_valence() ;   

    // Les verifications :
    assert ((ind_1 >= 0) && (ind_1 < val)) ;
    assert ((ind_2 >= 0) && (ind_2 < val)) ;
    assert (source.get_index_type(ind_1) != source.get_index_type(ind_2)) ;

    // On veut ind_1 < ind_2 :
    if (ind_1 > ind_2) {
		int auxi = ind_2 ;
		ind_2 = ind_1 ;
		ind_1 = auxi ;
    }
    
    // On construit le resultat :
    int val_res = val - 2 ;
   
    Itbl tipe(val_res) ;
	
    for (int i=0 ; i<ind_1 ; i++)
		tipe.set(i) = source.get_index_type(i) ;
    for (int i=ind_1 ; i<ind_2-1 ; i++)
		tipe.set(i) = source.get_index_type(i+1) ;
    for (int i = ind_2-1 ; i<val_res ; i++)
		tipe.set(i) = source.get_index_type(i+2) ;
	
    Tensor res(source.get_mp(), val_res, tipe, source.get_triad()) ; 
	
    Scalar work(source.get_mp()) ;
    
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_source(val) ;
	
    for (int i=0 ; i<res.get_n_comp() ; i++) {

		Itbl jeux_indice_res(res.indices(i)) ;
		
		for (int j=0 ; j<ind_1 ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j) ;
		for (int j=ind_1+1 ; j<ind_2 ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j-1) ;
		for (int j=ind_2+1 ; j<val ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j-2) ;
	    
		work.set_etat_zero() ;
		for (int j=1 ; j<=3 ; j++) {
	    	jeux_indice_source.set(ind_1) = j ;
	    	jeux_indice_source.set(ind_2) = j ;
	    	work = work + source(jeux_indice_source) ;
	    }
	    
		res.set(jeux_indice_res) = work ;
	}
	
    return res ;
}




				//------------------//
				//     diffrel	    //
				//------------------//


Tbl diffrel(const Tensor& aa, const Tensor& bb, ostream& ost) {

	int val = aa.get_valence() ; 

	assert(bb.get_valence() == val) ; 
	
	int n_comp_a = aa.get_n_comp() ; 
	int n_comp_b = bb.get_n_comp() ; 
	
	const Tensor* tmax ; 
	int n_comp_max ; 
	if (n_comp_a >= n_comp_b) {
		n_comp_max = n_comp_a ; 
		tmax = &aa ; 
	}
	else {
		n_comp_max = n_comp_b ; 
		tmax = &bb ; 
	}
	
	int nz = aa.get_mp().get_mg()->get_nzone() ; 
	Tbl resu(n_comp_max, nz) ; 
	resu.set_etat_qcq() ; 

	Itbl idx(val) ; 
	
	for (int ic=0; ic<n_comp_max; ic++) {
		idx = tmax->indices(ic) ; 
		Tbl diff = diffrel( aa(idx), bb(idx) ) ; 
		
		if (n_comp_max > 1) ost << "   Comp." ; 
		for (int j=0 ; j<val ; j++) {
	  		ost << " " << idx(j) ;
      	}
		if (n_comp_max > 1) ost << " : " ; 
                else ost << "   " ;
		for (int l=0; l<nz; l++) {
			ost << "  " << diff(l) ;
			resu.set(ic, l) = diff(l) ; 
		}
		ost << "\n" ; 
		
	}
	
	return resu ; 
}


				//--------------------//
				//     diffrelmax	  //
				//--------------------//


Tbl diffrelmax(const Tensor& aa, const Tensor& bb, ostream& ost) {

	int val = aa.get_valence() ; 

	assert(bb.get_valence() == val) ; 
	
	int n_comp_a = aa.get_n_comp() ; 
	int n_comp_b = bb.get_n_comp() ; 
	
	const Tensor* tmax ; 
	int n_comp_max ; 
	if (n_comp_a >= n_comp_b) {
		n_comp_max = n_comp_a ; 
		tmax = &aa ; 
	}
	else {
		n_comp_max = n_comp_b ; 
		tmax = &bb ; 
	}
	
	int nz = aa.get_mp().get_mg()->get_nzone() ; 
	Tbl resu(n_comp_max, nz) ; 
	resu.set_etat_qcq() ; 

	Itbl idx(val) ; 
	
	for (int ic=0; ic<n_comp_max; ic++) {
		idx = tmax->indices(ic) ; 
		Tbl diff = diffrelmax( aa(idx), bb(idx) ) ; 
		
		if (n_comp_max > 1) ost << "   Comp." ; 
		for (int j=0 ; j<val ; j++) {
	  		ost << " " << idx(j) ;
      	}
		if (n_comp_max > 1) ost << " : " ; 
                else ost << "   " ;
		for (int l=0; l<nz; l++) {
			ost << "  " << diff(l) ;
			resu.set(ic, l) = diff(l) ; 
		}
		ost << "\n" ; 
		
	}
	
	return resu ; 
}



				//----------------//
				//     max	  //
				//----------------//


Tbl max(const Tensor& aa, ostream& ost) {

	int val = aa.get_valence() ; 

	int n_comp = aa.get_n_comp() ; 
	
	int nz = aa.get_mp().get_mg()->get_nzone() ; 
	Tbl resu(n_comp, nz) ; 
	resu.set_etat_qcq() ; 

	Itbl idx(val) ; 
	
	for (int ic=0; ic<n_comp; ic++) {

		idx = aa.indices(ic) ; 
		Tbl diff = max( aa(idx) ) ; 
		
		if (val > 0) ost << "   Comp." ; 
		for (int j=0 ; j<val ; j++) {
	  		ost << " " << idx(j) ;
      	}
		if (val > 0) ost << " : " ; 
                else ost << "   " ;
		for (int l=0; l<nz; l++) {
			ost << "  " << diff(l) ;
			resu.set(ic, l) = diff(l) ; 
		}
		ost << "\n" ; 
		
	}
	
	return resu ; 
}



				//----------------//
				//     min	  //
				//----------------//


Tbl min(const Tensor& aa, ostream& ost) {

	int val = aa.get_valence() ; 

	int n_comp = aa.get_n_comp() ; 
	
	int nz = aa.get_mp().get_mg()->get_nzone() ; 
	Tbl resu(n_comp, nz) ; 
	resu.set_etat_qcq() ; 

	Itbl idx(val) ; 
	
	for (int ic=0; ic<n_comp; ic++) {

		idx = aa.indices(ic) ; 
		Tbl diff = min( aa(idx) ) ; 
		
		if (val > 0) ost << "   Comp." ; 
		for (int j=0 ; j<val ; j++) {
	  		ost << " " << idx(j) ;
      	}
		if (val > 0) ost << " : " ; 
                else ost << "   " ;
		for (int l=0; l<nz; l++) {
			ost << "  " << diff(l) ;
			resu.set(ic, l) = diff(l) ; 
		}
		ost << "\n" ; 
		
	}
	
	return resu ; 
}


				//--------------------//
				//     maxabs	      //
				//--------------------//


Tbl maxabs(const Tensor& aa, ostream& ost) {

	int val = aa.get_valence() ; 

	int n_comp = aa.get_n_comp() ; 
	
	int nz = aa.get_mp().get_mg()->get_nzone() ; 
	Tbl resu(n_comp, nz) ; 
	resu.set_etat_qcq() ; 

	Itbl idx(val) ; 
	
	for (int ic=0; ic<n_comp; ic++) {

		idx = aa.indices(ic) ; 
		Tbl diff = max( abs( aa(idx) ) ) ; 
		
		if (val > 0) ost << "   Comp." ; 
		for (int j=0 ; j<val ; j++) {
	  		ost << " " << idx(j) ;
      	        }
		if (val > 0 ) ost << " : " ; 
                else ost << "   " ; 
		for (int l=0; l<nz; l++) {
			ost << "  " << diff(l) ;
			resu.set(ic, l) = diff(l) ; 
		}
		ost << "\n" ; 
		
	}
	
	return resu ; 
}






