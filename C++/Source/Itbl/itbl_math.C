/*
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


char itbl_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.3  1999/11/25  13:02:12  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/11/25  12:42:12  phil
 * conversion double->int
 *
 * Revision 2.1  1999/11/24  09:32:01  eric
 * fabs -> abs
 *
 * Revision 2.0  1999/11/17  16:04:44  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
/*
 * Surcharge des
 * Fonctions mathematiques applicables aux classes de type
 * 
 *	Itbl
 *
 * Typiquement: max, norme ...
 *
 */
 
// Headers C
// ---------
#include <math.h>
#include <stdlib.h>

// Headers Lorene
// --------------
#include "itbl.h"


			    //----------------//
			    // Absolute value //
			    //----------------//

Itbl abs(const Itbl& ti)
{
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas ETATZERO
    if (ti.get_etat() == ETATZERO) {
	return ti ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;	// sinon...

    Itbl to(ti.dim) ;		    // Itbl resultat
    to.set_etat_qcq() ;

    const int* xi = ti.t ; 
    int* xo = to.t ; 
    int taille = ti.get_taille() ;

    for (int i=0 ; i<taille ; i++) {
	xo[i] = abs( xi[i] ) ;
    }
    
    return to ;
}

		    //-------------------------------//
    	    	    //            max                //
		    //-------------------------------//

int max(const Itbl& ti) {
    
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (ti.get_etat() == ETATZERO) {
    	return 0 ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;     // sinon....
    
    const int* x = ti.t ; 
    int resu = x[0] ;
    for (int i=1; i<ti.get_taille(); i++) {
    	if ( x[i] > resu ) resu = x[i] ;
    }
	
    return resu ; 
}

		    //-------------------------------//
    	    	    //            min                //
		    //-------------------------------//

int min(const Itbl& ti) {
    
    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    // Cas particulier
    if (ti.get_etat() == ETATZERO) {
    	return 0 ;
    }
    
    // Cas general
    assert(ti.get_etat() == ETATQCQ) ;     // sinon....
    
    const int* x = ti.t ; 
    int resu = x[0] ;
    for (int i=1; i<ti.get_taille(); i++) {
    	if ( x[i] < resu ) resu = x[i] ;
    }
	
    return resu ; 
}

		    //-------------------------------//
    	    	    //            norme              //
		    //-------------------------------//

int norme(const Itbl& ti) {

    // Protection
    assert(ti.get_etat() != ETATNONDEF) ;
    
    int resu = 0 ; 
    
    if (ti.get_etat() != ETATZERO) {  // on n'effectue la somme que si necessaire

	assert(ti.get_etat() == ETATQCQ) ;     // sinon....
	const int* x = ti.t ; 
	for (int i=0; i<ti.get_taille(); i++) {
	    resu += abs( x[i] ) ;
	}

    }

    return resu ;
}

		    //-------------------------------//
    	    	    //          diffrel              //
		    //-------------------------------//

double diffrel(const Itbl& t1, const Itbl& t2) {
    
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    
    int norm2 = norme(t2) ;
    int normdiff = norme(t1-t2) ;
    double resu ;  
    if ( norm2 == 0 ) {
	resu = double(normdiff) ; 
    }
    else {
	resu =  double(normdiff) / double(norm2) ; 
    }
    
    return resu ; 
    
}

		    //-------------------------------//
    	    	    //       diffrelmax              //
		    //-------------------------------//

double diffrelmax(const Itbl& t1, const Itbl& t2) {
    
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    
    int max2 = max(abs(t2)) ;
    int maxdiff = max(abs(t1-t2)) ;
    double resu ;  
    if ( max2 == 0 ) {
	resu = double(maxdiff) ; 
    }
    else {
	resu =  double(maxdiff) / double(max2) ; 
    }
    
    return resu ; 
    
}
