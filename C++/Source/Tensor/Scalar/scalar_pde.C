/*
 * Methods of the class Scalar for various partial differential equations
 *
 *  See file scalar.h for documentation. 
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Eric Gourgoulhon (for preceding class Cmp)
 *   Copyright (c) 1999-2001 Philippe Grandclement (for preceding class Cmp)
 *   Copyright (c) 2000-2001 Jerome Novak (for preceding class Cmp)
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


char scalar_pde_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/25 08:06:56  e_gourgoulhon
 * First versions (use Cmp as intermediate quantities).
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "tensor.h"
#include "param.h"
#include "cmp.h"

		    //-----------------------------------//
		    //      Scalar Poisson equation	 //
		    //-----------------------------------//

// Version without parameters
// --------------------------

Scalar Scalar::poisson() const {
    
    Param bidon ;
	Cmp csource(*this) ; 
    Cmp cresu(mp) ;     
	
    mp->poisson(csource, bidon, cresu) ; 

	Scalar resu(cresu) ; 
    return resu ;          
}

// Version with parameters
// -----------------------

void Scalar::poisson(Param& par, Scalar& uu) const {
    
	Cmp csource(*this) ; 
    Cmp cuu(mp) ;     

    mp->poisson(csource, par, cuu) ;     
    
	uu = cuu ; 
}



		    //-----------------------------------//
		    //      Scalar d'Alembert equation	 //
		    //-----------------------------------//

Scalar Scalar::avance_dalembert(Param& par, const Scalar& fjm1, 
	const Scalar& source) const {
  
	Cmp csource(source) ; 
	Cmp cfjm1(fjm1) ; 
	Cmp cfjp1(mp) ;
	Cmp cuu(*this) ; 
  
  	mp->dalembert(par, cfjp1, cuu, cfjm1, csource) ;

	Scalar fjp1(cfjp1) ; 
	
	return fjp1 ;
}
