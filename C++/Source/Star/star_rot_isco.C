/*
 * Method of class Star_rot to compute the location of the ISCO
 *
 * (see file star_rot.h for documentation).
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 J. Leszek Zdunik
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


char star_rot_isco_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "star_rot.h"
#include "param.h"
#include "utilitaires.h"

double fonct_etoile_rot_isco(double, const Param&) ;


//=============================================================================
//		r_isco()
//=============================================================================

double Star_rot::r_isco(ostream* ost) const {

    if (p_r_isco == 0x0) {    // a new computation is required

    cout << "Star_rot::r_isco : not implemented yet !" << endl ; 
    abort() ; 

    }  // End of computation

    return *p_r_isco ;

}



//=============================================================================
//		f_isco()
//=============================================================================

double Star_rot::f_isco() const {

    if (p_f_isco == 0x0) {    // a new computation is required

    	r_isco() ; 		// f_isco is computed in the method r_isco()

    	assert(p_f_isco != 0x0) ;
    }

    return *p_f_isco ;

}

//=============================================================================
//		lspec_isco()
//=============================================================================

double Star_rot::lspec_isco() const {

    if (p_lspec_isco == 0x0) {    // a new computation is required

    	r_isco() ; 	// lspec_isco is computed in the method r_isco()

    	assert(p_lspec_isco != 0x0) ;
    }

    return *p_lspec_isco ;

}

//=============================================================================
//		espec_isco()
//=============================================================================

double Star_rot::espec_isco() const {

    if (p_espec_isco == 0x0) {    // a new computation is required

    	r_isco() ; 	// espec_isco is computed in the method r_isco()

    	assert(p_espec_isco != 0x0) ;
    }

    return *p_espec_isco ;

}


//=============================================================================
//              f_eq()
//=============================================================================

double Star_rot::f_eq() const {
  
  if (p_f_eq == 0x0) {    // a new computation is required

    r_isco() ;              // f_eq is computed in the method r_isco()

    assert(p_f_eq != 0x0) ;
  }

  return *p_f_eq ;

}







