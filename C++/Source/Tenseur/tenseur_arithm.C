/*
 *  Arithmetics functions for the Tenseur class.
 *
 *  These functions are not member functions of the Tenseur class.
 *
 *  (see file tenseur.h for documentation).
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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

char tenseur_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2014/10/13 08:53:42  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:13:18  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2003/10/13 10:33:43  f_limousin
 * *** empty log message ***
 *
 * Revision 1.6  2003/06/20 14:52:21  f_limousin
 * Put an assert on "poids" into comments
 *
 * Revision 1.5  2003/03/03 19:40:52  f_limousin
 * Suppression of an  assert on a metric associated with a tensor.
 *
 * Revision 1.4  2002/10/16 14:37:14  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/09/06 14:49:25  j_novak
 * Added method lie_derive for Tenseur and Tenseur_sym.
 * Corrected various errors for derive_cov and arithmetic.
 *
 * Revision 1.2  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  2000/02/09  19:30:22  eric
 * MODIF IMPORTANTE: la triade de decomposition est desormais passee en
 * argument des constructeurs.
 *
 * Revision 2.4  2000/02/08  19:04:58  eric
 * Les fonctions arithmetiques ne sont plus amies.
 * Les fonctions exp, log et sqrt se trouvent desormais dans le fichier
 *   tenseur_math.C
 * Modif de diverses operations (notament division avec double)
 * Ajout de nouvelles operations (par ex. Tenseur + double, etc...)
 *
 * Revision 2.3  2000/02/01  15:40:29  eric
 * Ajout de la fonction sqrt
 *
 * Revision 2.2  2000/01/11  11:14:15  eric
 * Changement de nom pour la base vectorielle : base --> triad
 *
 * Revision 2.1  2000/01/10  17:25:34  eric
 * Gestion des bases vectorielles (triades de decomposition).
 *
 * Revision 2.0  1999/12/02  17:18:47  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>
#include <cmath>

// Headers Lorene
#include "tenseur.h"

			//********************//
			// OPERATEURS UNAIRES //
			//********************//

namespace Lorene {
Tenseur operator+(const Tenseur & t) {

    return t ; 

}

Tenseur operator-(const Tenseur & t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    if (t.get_etat() == ETATZERO)
	return t ;
    else { 
	Tenseur res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
		    t.get_triad(), t.get_metric(), t.get_poids()) ;


	res.set_etat_qcq();

	for (int i=0 ; i<res.get_n_comp() ; i++) {
	    Itbl indices (res.donne_indices(i)) ;    
	    res.set(indices) = -t(indices) ;
	    }
	return res ;
	}
}

			//**********//
			// ADDITION //
			//**********//

Tenseur operator+(const Tenseur & t1, const Tenseur & t2) {
    
    assert ((t1.get_etat() != ETATNONDEF) && (t2.get_etat() != ETATNONDEF)) ;
    assert (t1.get_valence() == t2.get_valence()) ;
    assert (t1.get_mp() == t2.get_mp()) ;
    if (t1.get_valence() != 0) {
	assert ( *(t1.get_triad()) == *(t2.get_triad()) ) ;
    }
    
    for (int i=0 ; i<t1.get_valence() ; i++)
	assert(t1.get_type_indice(i) == t2.get_type_indice(i)) ;
    //    assert (t1.get_metric() == t2.get_metric()) ;
    //assert (fabs(t1.get_poids() - t2.get_poids())<1.e-10) ;

    if (t1.get_etat() == ETATZERO)
	return t2 ;
    else if (t2.get_etat() == ETATZERO)
	    return t1 ;
	 else {
	    Tenseur res(*(t1.get_mp()), t1.get_valence(), t1.get_type_indice(), 
			t1.get_triad(), t1.get_metric(), t1.get_poids() ) ;

	    res.set_etat_qcq() ;
	    for (int i=0 ; i<res.get_n_comp() ; i++) {
		Itbl indices (res.donne_indices(i)) ;
		res.set(indices) = t1(indices) + t2(indices) ;
		}
	return res ;
	}
}


Tenseur operator+(const Tenseur & t1, double x) {
    
    assert (t1.get_etat() != ETATNONDEF) ;
    assert (t1.get_valence() == 0) ;
    
    if (x == double(0)) {
	return t1 ;
    }
    
    Tenseur res( *(t1.get_mp()), t1.get_metric(), t1.get_poids() ) ;

    res.set_etat_qcq() ;

    res.set() = t1() + x ;	// Cmp + double

    return res ;

}


Tenseur operator+(double x, const Tenseur & t2) {
    
    return t2 + x ; 
    
}

Tenseur operator+(const Tenseur & t1, int m) {
    
    return t1 + double(m) ;

}


Tenseur operator+(int m, const Tenseur & t2) {
    
    return t2 + double(m) ; 
    
}



			//**************//
			// SOUSTRACTION //
			//**************//

Tenseur operator-(const Tenseur & t1, const Tenseur & t2) {

    return (t1 + (-t2)) ;

}


Tenseur operator-(const Tenseur & t1, double x) {

    assert (t1.get_etat() != ETATNONDEF) ;
    assert (t1.get_valence() == 0) ;
    
    if (x == double(0)) {
	return t1 ;
    }
    
    Tenseur res( *(t1.get_mp()), t1.get_metric(), t1.get_poids() ) ;

    res.set_etat_qcq() ;

    res.set() = t1() - x ;	// Cmp - double

    return res ;

}


Tenseur operator-(double x, const Tenseur & t2) {
    
    return - (t2 - x) ; 
    
}


Tenseur operator-(const Tenseur & t1, int m) {

    return t1 - double(m) ; 
    
}


Tenseur operator-(int m, const Tenseur & t2) {

    return - (t2 - double(m)) ;     

}



			//****************//
			// MULTIPLICATION //
			//****************//



Tenseur operator*(double x, const Tenseur& t) {
    
    assert (t.get_etat() != ETATNONDEF) ;
    if ( (t.get_etat() == ETATZERO) || (x == double(1)) )
	return t ;
    else {
	Tenseur res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
		    t.get_triad(), t.get_metric(), t.get_poids() ) ;

	if ( x == double(0) )
	    res.set_etat_zero() ;
	else {
	    res.set_etat_qcq() ;
	    for (int i=0 ; i<res.get_n_comp() ; i++) {
		Itbl indices (res.donne_indices(i)) ;
		res.set(indices) = x*t(indices) ;
		}
	    }
	    return res ; 
	}
}


Tenseur operator* (const Tenseur& t, double x) {
    return x * t ;
}

Tenseur operator*(int m, const Tenseur& t) {
    return double(m) * t ; 
}


Tenseur operator* (const Tenseur& t, int m) {
    return double(m) * t ;
}


			//**********//
			// DIVISION //
			//**********//

Tenseur operator/ (const Tenseur& t1, const Tenseur& t2) {
    
    // Protections
    assert(t1.get_etat() != ETATNONDEF) ;
    assert(t2.get_etat() != ETATNONDEF) ;
    assert(t2.get_valence() == 0) ; // t2 doit etre un scalaire !
    assert(t1.get_mp() == t2.get_mp()) ;

    double poids_res = t1.get_poids() - t2.get_poids() ;
    poids_res = (fabs(poids_res) < 1.e-10 ? 0. : poids_res) ;
    const Metrique* met_res = 0x0 ;
    if (poids_res != 0.) {
      //    assert((t1.get_metric() != 0x0) || (t2.get_metric() != 0x0)) ;
      if (t1.get_metric() != 0x0) met_res = t1.get_metric() ;
      else met_res = t2.get_metric() ;
    }
    // Cas particuliers
    if (t2.get_etat() == ETATZERO) {
	cout << "Division by 0 in Tenseur / Tenseur !" << endl ;
	abort() ; 
    }
    if (t1.get_etat() == ETATZERO) {
        Tenseur resu(t1) ;
	resu.set_poids(poids_res) ;
	resu.set_metric(*met_res) ;
    	return resu ;
    }

    // Cas general
    
    assert(t1.get_etat() == ETATQCQ) ;  // sinon...
    assert(t2.get_etat() == ETATQCQ) ;  // sinon...

    Tenseur res(*(t1.get_mp()), t1.get_valence(), t1.get_type_indice(), 
		t1.get_triad(), met_res, poids_res) ;

    res.set_etat_qcq() ;
    for (int i=0 ; i<res.get_n_comp() ; i++) {
	Itbl indices (res.donne_indices(i)) ;
	res.set(indices) = t1(indices) / t2() ;	    // Cmp / Cmp
    }
    return res ;

}


Tenseur operator/ (const Tenseur& t, double x) {

    assert (t.get_etat() != ETATNONDEF) ;
 
    if ( x == double(0) ) {
	cout << "Division by 0 in Tenseur / double !" << endl ;
	abort() ;
    }

    if ( (t.get_etat() == ETATZERO) || (x == double(1)) )
	return t ;
    else {
	Tenseur res(*(t.get_mp()), t.get_valence(), t.get_type_indice(), 
		    t.get_triad(), t.get_metric(), t.get_poids()) ;

	res.set_etat_qcq() ;
	for (int i=0 ; i<res.get_n_comp() ; i++) {
	    Itbl indices (res.donne_indices(i)) ;
	    res.set(indices) = t(indices) / x ;	    // Cmp / double
	}
	return res ; 
    }

}




Tenseur operator/ (double x, const Tenseur& t) {
    
    if (t.get_etat() == ETATZERO) {
	cout << "Division by 0 in double / Tenseur !" << endl ;
	abort() ; 
    }
    
    assert (t.get_etat() == ETATQCQ) ;
    assert(t.get_valence() == 0) ;	// Utilisable que sur scalaire !
    
    Tenseur res( *(t.get_mp()), t.get_metric(), -t.get_poids() ) ;
    res.set_etat_qcq() ;
    res.set() = x / t() ;	// double / Cmp
    return res ;
}


Tenseur operator/ (const Tenseur& t, int m) {

    return t / double(m) ; 
}


Tenseur operator/ (int m, const Tenseur& t) {
    
    return double(m) / t ; 
}





}
