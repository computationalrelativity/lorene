/*
 *  Methods of the class Cmp for partial differential equations
 *   with a falloff condition at the outer boundary
 *
 *    (see file cmp.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Joshua A. Faber
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

char cmp_pde_falloff_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/11/30 20:47:38  k_taniguchi
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"
#include "param.h"

		    //-----------------------------------//
		    //      Scalar Poisson equation	 //
		    //-----------------------------------//

// Version without parameters
// --------------------------

Cmp Cmp::poisson_falloff(int k_falloff) const {
    
    Param bidon ;
    Cmp resu(*mp) ; 
    
    mp->poisson_falloff(*this, bidon, resu, k_falloff) ; 

    return resu ;          
}

// Version with parameters
// -----------------------

void Cmp::poisson_falloff(Param& par, Cmp& uu, int k_falloff) const {
    
    mp->poisson_falloff(*this, par, uu, k_falloff) ;     
    
}
