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
 * Revision 1.3  2002/08/08 15:10:45  j_novak
 * The flag "plat" has been added to the class Metrique to show flat metrics.
 *
 * Revision 1.2  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
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

    Tenseur res( *(t.get_mp()), t.get_metric(), 0.5*t.get_poids() ) ;
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
    
    Tenseur res( *(t.get_mp()), t.get_metric(), t.get_poids() ) ;//!?
    res.set_etat_qcq() ;
    res.set() = abs(t()) ;
    return res ;
}

			    //--------------//
			    // Power^double //
			    //--------------//

Tenseur pow (const Tenseur& t, double a) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
        
    Tenseur res( *(t.get_mp()), t.get_metric(), a*t.get_poids() ) ;
    res.set_etat_qcq() ;
    if (t.get_etat() == ETATZERO)
      if (a > double(0)) {
	res.set_etat_zero() ;
	res.set_std_base() ;
	return res ; 
      }
      else {
	cout << "pow(Tenseur, double) : ETATZERO^x with x <= 0 !" << endl ; 
	abort() ;
      } 
    else {
      assert(t.get_etat() == ETATQCQ) ;
      res.set() = pow( t(), a ) ;
    }
    res.set_std_base() ;
    return res ;
}

			    //--------------//
			    //   Power^int  //
			    //--------------//

Tenseur pow (const Tenseur& t, int n) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    assert (t.get_valence() == 0) ;	// Scalaire uniquement ...
        
    Tenseur res( *(t.get_mp()), t.get_metric(), n*t.get_poids() ) ;
    res.set_etat_qcq() ;
    if (t.get_etat() == ETATZERO)
      if (n > double(0)) {
	res.set_etat_zero() ;
	return res ; 
      }
      else {
	cout << "pow(Tenseur, int) : ETATZERO^n with n <= 0 !" << endl ; 
	abort() ;
      } 
    else {
      assert(t.get_etat() == ETATQCQ) ;
      res.set() = pow( t(), n ) ;
    }
    return res ;
}

