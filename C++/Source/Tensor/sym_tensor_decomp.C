/*
 *  Methods transverse( ) and longit_pot( ) of class Sym_tensor
 *
 *    (see file sym_tensor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon & Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char sym_tensor_decomp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/11/27 16:01:47  e_gourgoulhon
 * First implmentation.
 *
 * Revision 1.1  2003/11/26 21:57:03  e_gourgoulhon
 * First version; not ready yet.
 *
 *
 * $Header$
 *
 */


// C headers
#include <stdlib.h>

// Lorene headers
#include "tensor.h"
#include "metric.h"

const Sym_tensor_trans& Sym_tensor::transverse(const Metric& metre) const {

	set_dependance(metre) ;
	int jp = get_place_met(metre) ;
	assert ((jp>=0) && (jp<N_MET_MAX)) ;

	if (p_transverse[jp] == 0x0) { // a new computation is necessary

		assert( (type_indice(0) == CON) && (type_indice(1) == CON) ) ; 

		const Vector& ww = longit_pot(metre) ;
		
		const Tensor& dww = ww.derive_con(metre) ; 
		
		p_transverse[jp] = new Sym_tensor_trans(*mp, *triad, metre) ;
		
		for (int i=1; i<=3; i++) {
			for (int j=i; j<=3; j++) {
				p_transverse[jp]->set(i,j) = operator()(i,j) 
					- dww(i,j) - dww(j,i) ; 
			}
		}

	}

  	return *p_transverse[jp] ;
	

}




const Vector& Sym_tensor::longit_pot(const Metric& metre) const {

	set_dependance(metre) ;
	int jp = get_place_met(metre) ;
	assert ((jp>=0) && (jp<N_MET_MAX)) ;

	if (p_longit_pot[jp] == 0x0) {  // a new computation is necessary
		
		const Metric_flat* metf = dynamic_cast<const Metric_flat*>(&metre) ; 
		if (metf == 0x0) {
			cout << "Sym_tensor::longit_pot : the case of a non flat metric"
			 << endl <<"  is not treated yet !" << endl ; 
			abort() ; 
		}
		
		Vector hhh = divergence(metre) ; 
		for (int i=1; i<=3; i++) {
			hhh.set(i).inc_dzpuis(2) ; 
		}
				
		p_longit_pot[jp] = new Vector( hhh.poisson(double(1)) ) ; 
		
	}

	return *p_longit_pot[jp] ;
	

}

