/*
 * Computations of partial derivatives d/dx, d/dy and d/dz of a Cmp.
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char cmp_deriv_C[] = "$Header$" ;



/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.7  2000/09/11  15:55:12  eric
 * Calcul de dsdx, dsdy et dsdz: suppression des methodes Map::deriv_x, etc...
 *  et introduction des methodes Map::comp_x_from_spherical, etc...
 *
 * Revision 1.6  2000/02/08  14:20:14  phil
 * vire set_etat_qcq
 *
 * Revision 1.5  2000/01/27  17:27:49  phil
 * coorection etat du resultat
 *
 * Revision 1.4  2000/01/26  13:11:36  eric
 * Modifs pour tenir compte du reprototypage complet des routines de derivation
 * des mappings (Map::dsdr, etc...). Le resultat p_* est desormais alloue
 * a l'exterieur de la routine Map::*.
 *
 * Revision 1.3  1999/11/29  14:38:16  eric
 * *** empty log message ***
 *
 * Revision 1.2  1999/11/29  12:57:07  eric
 * Introduction du laplacien.
 *
 * Revision 1.1  1999/11/25  16:28:17  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */
 
// Headers C
#include <cstdlib>

// Headers Lorene
#include "cmp.h"

			//---------------------//
			//	d/dr	       //
			//---------------------//

namespace Lorene {
const Cmp& Cmp::dsdr() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdr == 0x0) {
	p_dsdr = new Cmp(mp) ; 
        mp->dsdr(*this, *p_dsdr) ;
    }
    
    return *p_dsdr ;

}

			//------------------------//
			//	1/r d/dtheta      //
			//------------------------//

const Cmp& Cmp::srdsdt() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srdsdt == 0x0) {
	p_srdsdt = new Cmp(mp) ;
        mp->srdsdt(*this, *p_srdsdt) ;
    }
    
    return *p_srdsdt ;

}


			//----------------------------------//
			//	1/(r sin(theta) d/dphi	    //
			//----------------------------------//

const Cmp& Cmp::srstdsdp() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srstdsdp == 0x0) {
        p_srstdsdp = new Cmp(mp) ; 
	mp->srstdsdp(*this, *p_srstdsdp) ;
    }
    
    return *p_srstdsdp ;

}

			//---------------------//
			//	d/dx	       //
			//---------------------//

const Cmp& Cmp::dsdx() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdx == 0x0) {
        p_dsdx = new Cmp(mp) ; 
	mp->comp_x_from_spherical(dsdr(), srdsdt(), srstdsdp(), *p_dsdx) ;
    }
    
    return *p_dsdx ;

}

			//---------------------//
			//	d/dy	       //
			//---------------------//

const Cmp& Cmp::dsdy() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdy == 0x0) {
        p_dsdy = new Cmp(mp) ; 
	mp->comp_y_from_spherical(dsdr(), srdsdt(), srstdsdp(), *p_dsdy) ;
    }
    
    return *p_dsdy ;

}

			//---------------------//
			//	d/dz	       //
			//---------------------//

const Cmp& Cmp::dsdz() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdz == 0x0) {
        p_dsdz = new Cmp(mp) ; 
	mp->comp_z_from_spherical(dsdr(), srdsdt(), *p_dsdz) ;
    }
    
    return *p_dsdz ;

}

			//---------------------//
			//	d/dx^i	       //
			//---------------------//

const Cmp& Cmp::deriv(int i) const {
    
    switch (i) {
	
	case 0 : {
	    return dsdx() ; 
	}
	
	case 1 : {
	    return dsdy() ; 
	}
	
	case 2 : {
	    return dsdz() ; 
	}
	
	default : {
	    cout << "Cmp::deriv : index i out of range !" << endl ; 
	    cout << "  i = " << i << endl ; 
	    abort() ; 
	    return dsdx() ;  // Pour satisfaire le compilateur !
	}
	
    }
    
}

			//---------------------//
			//     Laplacian       //
			//---------------------//

const Cmp& Cmp::laplacien(int zec_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the Laplacian has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 
    if ( (p_lap == 0x0) || (zec_mult_r != ind_lap) ) {
	if (p_lap != 0x0) {
	    delete p_lap ;  // the Laplacian had been computed but with
			    //  a different value of zec_mult_r
	}
	p_lap = new Cmp(mp) ;
	mp->laplacien(*this, zec_mult_r, *p_lap) ;
	ind_lap = zec_mult_r ;
    }
    
    return *p_lap ;
    
}
    
   
    
    
}
