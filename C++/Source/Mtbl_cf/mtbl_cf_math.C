/*
 * Mathematical functions related to the Mtbl_cf class
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char mtbl_cf_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:53:08  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:15  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.3  2000/02/25  10:57:40  eric
 * Suppressions des appels a nettoie().
 *
 * Revision 1.2  1999/10/29  15:46:20  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/10/29  15:08:13  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
// ---------
#include <cmath>
#include <cstdlib>

// Headers Lorene
// --------------
#include "mtbl_cf.h"


			    //----------------//
			    // Absolute value //
			    //----------------//

namespace Lorene {
Mtbl_cf abs(const Mtbl_cf& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general

    assert(ti.get_etat() == ETATQCQ) ;	// sinon...

    Mtbl_cf to(ti.get_mg(), ti.base) ;			// Mtbl_cf resultat
   
    to.set_etat_qcq() ;

    int nzone = ti.get_nzone() ;

    for (int l=0 ; l<nzone ; l++) {
	*(to.t[l]) = abs( *(ti.t[l]) ) ;	
    }

    return to ;
}




		    //-------------------------------//
    	    	    //            max                //
		    //-------------------------------//

Tbl max(const Mtbl_cf& mti) {

    // Protection
    assert(mti.get_etat() != ETATNONDEF) ;
    
    int nz = mti.get_nzone() ;
    
    Tbl resu(nz) ; 
    
    if (mti.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {  	// Cas general

	assert(mti.get_etat() == ETATQCQ) ;     // sinon....

	resu.set_etat_qcq() ; 
	for (int l=0 ; l<nz ; l++) {
	    resu.set(l) = max( *(mti.t[l]) ) ;
	}
    }
     
    return resu ; 
}

		    //-------------------------------//
    	    	    //            min                //
		    //-------------------------------//

Tbl min(const Mtbl_cf& mti) {

    // Protection
    assert(mti.get_etat() != ETATNONDEF) ;
    
    int nz = mti.get_nzone() ;
    
    Tbl resu(nz) ; 
    
    if (mti.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {  	// Cas general

	assert(mti.get_etat() == ETATQCQ) ;     // sinon....
    
	resu.set_etat_qcq() ; 
	for (int l=0 ; l<nz ; l++) {
	    resu.set(l) = min( *(mti.t[l]) ) ;
	}
    }
     
    return resu ; 
}

		    //-------------------------------//
    	    	    //            norme              //
		    //-------------------------------//

Tbl norme(const Mtbl_cf& mti) {

    // Protection
    assert(mti.get_etat() != ETATNONDEF) ;
    
    int nz = mti.get_nzone() ;
    
    Tbl resu(nz) ; 
    
    if (mti.get_etat() == ETATZERO) {
	resu.annule_hard() ; 
    }
    else {  	// Cas general

	assert(mti.get_etat() == ETATQCQ) ;     // sinon....
    
	resu.set_etat_qcq() ; 
	for (int l=0 ; l<nz ; l++) {
	    resu.set(l) = norme( *(mti.t[l]) ) ;
	}
    }
     
    return resu ; 
}

		    //-------------------------------//
    	    	    //          diffrel              //
		    //-------------------------------//

Tbl diffrel(const Mtbl_cf& mt1, const Mtbl_cf& mt2) {
    
    // Protections
    assert(mt1.get_etat() != ETATNONDEF) ;
    assert(mt2.get_etat() != ETATNONDEF) ;
    
    int nz = mt1.get_nzone() ;
    Tbl resu(nz) ; 
    
    Mtbl_cf diff = mt1 - mt2 ; 
    
    Tbl normdiff = norme(diff) ; 
    Tbl norme2 = norme(mt2) ;
    
    assert(normdiff.get_etat() == ETATQCQ) ;     
    assert(norme2.get_etat() == ETATQCQ) ; 

    resu.set_etat_qcq() ; 
    for (int l=0; l<nz; l++) {
	if ( norme2(l) == double(0) ) {
	    resu.set(l) = normdiff(l) ; 
	}
	else{
	    resu.set(l) = normdiff(l) / norme2(l) ; 		    
	}
    }
    
    return resu ; 
    
}

		    //-------------------------------//
    	    	    //          diffrelmax           //
		    //-------------------------------//

Tbl diffrelmax(const Mtbl_cf& mt1, const Mtbl_cf& mt2) {
    
    // Protections
    assert(mt1.get_etat() != ETATNONDEF) ;
    assert(mt2.get_etat() != ETATNONDEF) ;
    
    int nz = mt1.get_nzone() ;
    Tbl resu(nz) ; 
    
    Tbl max2 = max(abs(mt2)) ;

    Mtbl_cf diff = mt1 - mt2 ; 
    Tbl maxdiff = max(abs(diff)) ; 
    
    assert(maxdiff.get_etat() == ETATQCQ) ;     
    assert(max2.get_etat() == ETATQCQ) ; 

    resu.set_etat_qcq() ; 
    for (int l=0; l<nz; l++) {
	if ( max2(l) == double(0) ) {
	    resu.set(l) = maxdiff(l) ; 
	}
	else{
	    resu.set(l) = maxdiff(l) / max2(l) ; 		    
	}
    }
    
    return resu ; 
    
}
}
