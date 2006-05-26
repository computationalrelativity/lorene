/*
 *  Member functions of the class Map_radial for various theta manipulations
 *  of Scalar's.
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char map_radial_th_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2006/05/26 09:00:11  j_novak
 * New members for multiplication or division by cos(theta).
 *
 * Revision 1.2  2003/11/05 15:27:52  e_gourgoulhon
 * Treatment of case ETATUN now possible by a simple call
 *  to set_etat_qcq(), thanks to a modification of the latter.
 *
 * Revision 1.1  2003/11/04 22:59:13  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "map.h"
#include "tensor.h"




			//---------------------------//
			//          mult_cost        //
			//---------------------------//

void Map_radial::mult_cost(Scalar& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
		return ;			 // Nothing to do if the Scalar is null 
    }

    assert((ci.get_etat() == ETATQCQ) || (ci.get_etat() == ETATUN)) ;
	
    Valeur& val = ci.set_spectral_va() ; 

    assert(val.get_mg() == mg) ; 
     
	val = val.mult_ct() ; 	// Multiplication by cos(theta) 

	ci.set_etat_qcq() ; 
	
}

			//---------------------------//
			//          div_cost         //
			//---------------------------//

void Map_radial::div_cost(Scalar& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
		return ;			 // Nothing to do if the Scalar is null 
    }

    assert((ci.get_etat() == ETATQCQ) || (ci.get_etat() == ETATUN)) ;
            
    Valeur& val = ci.set_spectral_va() ; 

    assert(val.get_mg() == mg) ; 
         
    val = val.scost() ;		// Division by cos(theta)

	ci.set_etat_qcq() ; 
}

            
			//---------------------------//
			//          mult_sint        //
			//---------------------------//

void Map_radial::mult_sint(Scalar& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
		return ;			 // Nothing to do if the Scalar is null 
    }

    assert((ci.get_etat() == ETATQCQ) || (ci.get_etat() == ETATUN)) ;
	
    Valeur& val = ci.set_spectral_va() ; 

    assert(val.get_mg() == mg) ; 
     
	val = val.mult_st() ; 	// Multiplication by sin(theta) 

	ci.set_etat_qcq() ;             

}

			//---------------------------//
			//          div_sint         //
			//---------------------------//

void Map_radial::div_sint(Scalar& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
		return ;			 // Nothing to do if the Scalar is null 
    }

    assert((ci.get_etat() == ETATQCQ) || (ci.get_etat() == ETATUN)) ;
            
    Valeur& val = ci.set_spectral_va() ; 

    assert(val.get_mg() == mg) ; 
         
    val = val.ssint() ;		// Division by sin(theta)

	ci.set_etat_qcq() ; 
}



			//---------------------------//
			//          div_tant         //
			//---------------------------//

void Map_radial::div_tant(Scalar& ci) const {
    
    assert(ci.get_etat() != ETATNONDEF) ;
    
    if (ci.get_etat() == ETATZERO) {
		return ;			 // Nothing to do if the Scalar is null 
    }

    assert((ci.get_etat() == ETATQCQ) || (ci.get_etat() == ETATUN)) ;
            
    Valeur& val = ci.set_spectral_va() ; 

    assert(val.get_mg() == mg) ; 
     
	val = val.mult_ct() ; 	// Multiplication by cos(theta) 
    
    val = val.ssint() ;		// Division by sin(theta)

	ci.set_etat_qcq() ; 
}

