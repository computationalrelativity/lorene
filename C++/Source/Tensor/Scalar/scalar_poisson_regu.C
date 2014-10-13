/*
 * Method of regularization of the source of Poisson equation
 *
 * (see file scalar.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 2000-2001 Keisuke Taniguchi (for preceding Cmp version)
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


char scalar_poisson_regu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2003/10/11 14:46:00  e_gourgoulhon
 * Line 65-67: changed the name of local variable "triad" to "triad0"
 * in order not to shadow the class member triad.
 *
 * Revision 1.1  2003/09/25 08:56:28  e_gourgoulhon
 * First version (uses Cmp and Tenseur as intermediate quantities).
 *
 *
 * $Header$
 *
 */

// Header Lorene
#include "tensor.h"
#include "cmp.h"
#include "tenseur.h"
#include "param.h"

//******************************************************************

namespace Lorene {

void Scalar::poisson_regular(int k_div, int nzet, double unsgam1, Param& par,
			  Scalar& uu, Scalar& uu_regu, Scalar& uu_div,
			  Tensor& duu_div,
			  Scalar& source_regu, Scalar& source_div) const {
			  
	Cmp csource(*this) ; 
	Cmp cuu(uu) ;
	Cmp cuu_regu(uu_regu) ; 
	Cmp cuu_div(uu_div) ; 
	Cmp csource_regu(source_regu) ; 
	Cmp csource_div(source_div) ; 

	const Base_vect* triad0 = duu_div.get_triad() ; 
	
	Tenseur cduu_div(*mp, 1, COV, *triad0) ; 
	cduu_div.set_etat_qcq() ;
	Itbl ind(1) ;
	ind.set_etat_qcq() ; 
	for (int i=0; i<3; i++) {
		ind.set(0) = i+1 ;
		Cmp tmp( duu_div(ind) ) ;  
		cduu_div.set(i) = tmp ;  
	}
	
    mp->poisson_regular(csource, k_div, nzet, unsgam1, par,
			cuu, cuu_regu, cuu_div, cduu_div,
			csource_regu, csource_div) ;
	
	uu = cuu ; 
	uu_regu = uu ; 

	for (int i=1; i<=3; i++) {
		ind.set(0) = i ;
		duu_div.set(ind) = cduu_div(i-1) ;  
	}
	
	source_regu = csource_regu ; 
	source_div = csource_div ; 

}
}
