/*
 *  Mathematical functions for the Qtenseur class.
 *
 *  These functions are not member functions of the Qtenseur class.
 *
 *  (see file qtenseur.h for documentation).
 *
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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


char qtenseur_math_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/09/19 09:52:43  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "qtenseur.h"

			    //--------------//
			    // Exponential  //
			    //--------------//

Qtenseur exp (const Qtenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Qtenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    if (t.get_etat() == ETATZERO)
	res.set() = 1 ;
    else 
	res.set() = exp( t() ) ;
    return res ;
}


			    //---------------------//
			    // Neperian logarithm  //
			    //---------------------//

Qtenseur log (const Qtenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Qtenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    res.set() = log(t()) ;
    return res ;
}

			    //-------------//
			    // Square root //
			    //-------------//

Qtenseur sqrt(const Qtenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Qtenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    res.set() = sqrt(t()) ;
    return res ;
}

			    //----------------//
			    // Absolute value //
			    //----------------//

Qtenseur abs(const Qtenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Qtenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    res.set() = abs(t()) ;
    return res ;
}

