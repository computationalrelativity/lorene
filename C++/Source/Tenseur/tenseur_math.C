/*
 *  Mathematical functions for the Tenseur class.
 *
 *  These functions are not member functions of the Tenseur class.
 *
 *  (see file tenseur.h for documentation).
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2001 Jerome Novak
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


char tenseur_math_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:30  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.1  2001/06/18  13:56:25  novak
 * Ajout de la fonction abs() pour les scalaires
 *
 * Revision 2.0  2000/02/08 19:06:32  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tenseur.h"

			    //--------------//
			    // Exponential  //
			    //--------------//

Tenseur exp (const Tenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Tenseur res( *(t.get_mp()) ) ;
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

Tenseur log (const Tenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Tenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    res.set() = log(t()) ;
    return res ;
}

			    //-------------//
			    // Square root //
			    //-------------//

Tenseur sqrt(const Tenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Tenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    res.set() = sqrt(t()) ;
    return res ;
}

			    //----------------//
			    // Absolute value //
			    //----------------//

Tenseur abs(const Tenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
    
    Tenseur res( *(t.get_mp()) ) ;
    res.set_etat_qcq() ;
    res.set() = abs(t()) ;
    return res ;
}

