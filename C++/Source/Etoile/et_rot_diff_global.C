/*
 * Methods for computing global quantities within the class Et_rot_diff
 *
 * (see file et_rot_diff.h for documentation)
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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


char et_rot_diff_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/12/19 16:21:42  j_novak
 * Shadow hunt
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2001/10/24  16:06:53  eric
 * fonction tsw(): correction erreur ener. cin. (facteur 0.5).
 *
 * Revision 1.1  2001/10/19  08:18:30  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "et_rot_diff.h"
			//----------------------------//
			//	     T/W	      //
			//----------------------------//

double Et_rot_diff::tsw() const {

    if (p_tsw == 0x0) {    // a new computation is required
	
	Cmp dens = uuu() ; 

	dens.mult_r() ;			//  Multiplication by
	dens.va = (dens.va).mult_st() ;	//    r sin(theta)

	if (relativistic) {
	    dens = omega_field() * a_car() * b_car() * (ener_euler() + press()) 
			* dens ; 
	}
	else {    // Newtonian case 
	    dens = omega_field() * nbar() * dens ; 
	}

	dens.std_base_scal() ; 

	double tcin = 0.5 * dens.integrale() ;
	
	if (relativistic) {
	    
	    Cmp dens2 = a_car() * bbb() * gam_euler() * ener() ;
	    dens2.std_base_scal() ; 
	    double mass_p = dens2.integrale() ; 
	    
	    p_tsw = new double( tcin / ( mass_p + tcin - mass_g() ) ) ;  	
	   
	}
	else {	    // Newtonian case 
	    Cmp dens2 = 0.5 * nbar() * logn() ;
	    dens2.std_base_scal() ; 
	    double wgrav = dens2.integrale() ; 
	    p_tsw = new double( tcin / fabs(wgrav) ) ;  
	}


    }
    
    return *p_tsw ; 

}

