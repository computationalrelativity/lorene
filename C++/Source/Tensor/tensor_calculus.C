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
 * Revision 1.7  2004/02/18 18:50:07  e_gourgoulhon
 * -- Added new methods trace(...).
 * -- Tensorial product moved to file tensor_calculus_ext.C, since it is not
 *    a method of class Tensor.
 *
 * Revision 1.6  2004/02/18 15:54:23  e_gourgoulhon
 * Efficiency improved in method scontract: better handling of work (it is
 * now considered as a reference on the relevant component of the result).
 *
 * Revision 1.5  2003/12/05 16:38:50  f_limousin
 * Added method operator*
 *
 * Revision 1.4  2003/10/28 21:25:34  e_gourgoulhon
 * Method contract renamed scontract.
 *
 * Revision 1.3  2003/10/11 16:47:10  e_gourgoulhon
 * Suppressed the call to Ibtl::set_etat_qcq() after the construction
 * of the Itbl's, thanks to the new property of the Itbl class.
 *
 * Revision 1.2  2003/10/06 20:52:22  e_gourgoulhon
 * Added methods up, down and up_down.
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
#include "metric.h"


				//------------------//
				//       Trace      //
				//------------------//


Tensor Tensor::trace(int ind_1, int ind_2) const {
    
    // Les verifications :
    assert( (ind_1 >= 0) && (ind_1 < valence) ) ;
    assert( (ind_2 >= 0) && (ind_2 < valence) ) ;
    assert( ind_1 != ind_2 )  ;
    assert( type_indice(ind_1) != type_indice(ind_2) )  ;
 
    // On veut ind_1 < ind_2 :
    if (ind_1 > ind_2) {
		int auxi = ind_2 ;
		ind_2 = ind_1 ;
		ind_1 = auxi ;
    }
    
    // On construit le resultat :
    int val_res = valence - 2 ;
   
    Itbl tipe(val_res) ;
	
    for (int i=0 ; i<ind_1 ; i++)
		tipe.set(i) = type_indice(i) ;
    for (int i=ind_1 ; i<ind_2-1 ; i++)
		tipe.set(i) = type_indice(i+1) ;
    for (int i = ind_2-1 ; i<val_res ; i++)
		tipe.set(i) = type_indice(i+2) ;
	
    Tensor res(*mp, val_res, tipe, triad) ; 
	
    // Boucle sur les composantes de res :
	
    Itbl jeux_indice_source(valence) ;
	
    for (int i=0 ; i<res.get_n_comp() ; i++) {

		Itbl jeux_indice_res(res.indices(i)) ;
		
		for (int j=0 ; j<ind_1 ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j) ;
		for (int j=ind_1+1 ; j<ind_2 ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j-1) ;
		for (int j=ind_2+1 ; j<valence ; j++)
	    	jeux_indice_source.set(j) = jeux_indice_res(j-2) ;
	    
        Scalar& work = res.set(jeux_indice_res) ; 
		work.set_etat_zero() ;

		for (int j=1 ; j<=3 ; j++) {
	    	jeux_indice_source.set(ind_1) = j ;
	    	jeux_indice_source.set(ind_2) = j ;
	    	work += (*cmp[position(jeux_indice_source)]) ;
	    }
	    
	}
	
    return res ;
}


Tensor Tensor::trace(int ind1, int ind2, const Metric& gam) const {

    // Les verifications :
    assert( (ind1 >= 0) && (ind1 < valence) ) ;
    assert( (ind2 >= 0) && (ind2 < valence) ) ;
    assert( ind1 != ind2 )  ;
 
    if ( type_indice(ind1) != type_indice(ind2) ) {
        cout << "Tensor::trace(int, int, const Metric&) : Warning : \n"
            << "  the two indices for the trace have opposite types,\n"
            << "  hence the metric is useless !\n" ; 

        return trace(ind1, ind2) ;               
    }
    else {
        if ( type_indice(ind1) == COV  ) {
            return contract(gam.con(), 0, 1, *this, ind1, ind2) ; 
        }
        else{
            return contract(gam.cov(), 0, 1, *this, ind1, ind2) ; 
        }
    }
          
}



Scalar Tensor::trace() const {
    
    // Les verifications :
    assert( valence == 2 )  ;
    assert( type_indice(0) != type_indice(1) )  ;
 	
    Scalar res(*mp) ; 
	res.set_etat_zero() ; 
    
    for (int i=1; i<=3; i++) {
        res += operator()(i,i) ;    
    }
	
    return res ;
}


Scalar Tensor::trace(const Metric& gam) const {
    
    assert( valence == 2 )  ;
 	
    if ( type_indice(0) != type_indice(1) ) {
        cout << "Tensor::trace(const Metric&) : Warning : \n"
            << "  the two indices have opposite types,\n"
            << "  hence the metric is useless to get the trace !\n" ; 

        return trace() ;               
    }
    else {
        if ( type_indice(0) == COV  ) {
            return contract(gam.con(), 0, 1, *this, 0, 1) ; 
        }
        else{
            return contract(gam.cov(), 0, 1, *this, 0, 1) ; 
        }
    }
          
}


				//----------------------//
				//  Index manipulation	//
				//----------------------//


Tensor Tensor::up(int place, const Metric& met) const {
	
    assert (valence != 0) ;	    // Aucun interet pour un scalaire...
    assert ((place >=0) && (place < valence)) ;
    
    
	Tensor auxi = ::contract(met.con(), 1, *this, place) ;
    
    // On doit remettre les indices a la bonne place ...
    
    Itbl tipe(valence) ;

    for (int i=0 ; i<valence ; i++)
		tipe.set(i) = type_indice(i) ;
    tipe.set(place) = CON ;
    
    Tensor res(*mp, valence, tipe, triad) ;
    
    Itbl place_auxi(valence) ;
    
    for (int i=0 ; i<res.n_comp ; i++) {
	
		Itbl place_res(res.indices(i)) ;
	
		place_auxi.set(0) = place_res(place) ;
		for (int j=1 ; j<place+1 ; j++)
	    	place_auxi.set(j) = place_res(j-1)  ;
		place_res.set(place) = place_auxi(0) ;
	
		for (int j=place+1 ; j<valence ; j++)
			place_auxi.set(j) = place_res(j);	
	
		res.set(place_res) = auxi(place_auxi) ;
    }
	
    return res ;

} 


Tensor Tensor::down(int place, const Metric& met) const {
	
    assert (valence != 0) ;	    // Aucun interet pour un scalaire...
    assert ((place >=0) && (place < valence)) ;
    
	Tensor auxi = ::contract(met.cov(), 1, *this, place) ;
    
    // On doit remettre les indices a la bonne place ...
    
    Itbl tipe(valence) ;

    for (int i=0 ; i<valence ; i++)
		tipe.set(i) = type_indice(i) ;
    tipe.set(place) = COV ;
    
    Tensor res(*mp, valence, tipe, triad) ;
    
    Itbl place_auxi(valence) ;
    
    for (int i=0 ; i<res.n_comp ; i++) {
	
		Itbl place_res(res.indices(i)) ;
	
		place_auxi.set(0) = place_res(place) ;
		for (int j=1 ; j<place+1 ; j++)
	    	place_auxi.set(j) = place_res(j-1)  ;
		place_res.set(place) = place_auxi(0) ;
	
		for (int j=place+1 ; j<valence ; j++)
			place_auxi.set(j) = place_res(j);	
	
		res.set(place_res) = auxi(place_auxi) ;
    }
	
    return res ;

} 



Tensor Tensor::up_down(const Metric& met) const  {
    
    Tensor* auxi ;
    Tensor* auxi_old = new Tensor(*this) ;
    
    for (int i=0 ; i<valence ; i++) {

		if (type_indice(i) == COV) {
			auxi = new Tensor( auxi_old->up(i, met) ) ;
		}
		else{
			auxi = new Tensor( auxi_old->down(i, met) ) ;
		}
		
		delete auxi_old ;
		auxi_old = new Tensor(*auxi) ;
		delete auxi ;

    }
    
    Tensor result(*auxi_old) ;
    delete auxi_old ;

    return result ;
}














 
