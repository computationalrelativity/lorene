/*
 * Operator
 *
 *	f(x) |---->  (f(x)-f(1))/(x-1)
 *
 * in the external compactified domain,
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
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


 
 
/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:21  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:51  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/03/09  16:53:59  eric
 * Traitement du cas etat=ETATZERO
 *
 * Revision 2.2  1999/11/30  12:46:10  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.1  1999/10/18  13:43:19  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.0  1999/04/26  15:54:59  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cassert>

// Headers Lorene
#include "mtbl_cf.h"
#include "valeur.h"

// Local prototypes
namespace Lorene {
void _sxm1_identite(Tbl*, int&) ;
void _sxm1_cheb(Tbl *, int&) ;

// Version membre d'un Valeur
// --------------------------

void Valeur::sxm1_zec() {

    // Peut-etre ne rien faire ?
    if (etat==ETATZERO) {
	return ; 
    }

    assert(etat==ETATQCQ) ; 

    // Calcul des coef.
    coef() ;
    
    // Division par (x-1) dans la ZEC 
    c_cf->sxm1_zec() ;
    set_etat_cf_qcq() ;

    base = c_cf->base ; // On remonte la base de sortie au niveau Valeur
    
}


/*
 * Fonction membre de la classe Mtbl_cf pour la division par (x-1) 
 * dans la zone externe compactifiee applique a this
 * Par division par (x-1), il faut en fait entendre l'operateur
 * 
 *	f(x) |---->  (f(x)-f(1))/(x-1) 
 *
 *
 */

void Mtbl_cf::sxm1_zec()	   
{

// Routines de derivation
static void (*_sxm1[MAX_BASE])(Tbl *, int&) ;
static int nap = 0 ;

    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _sxm1[i] = _sxm1_identite ;
	}
	
	// Les routines existantes cas UNSURR
	_sxm1[R_CHEBU >> TRA_R] = _sxm1_cheb ;
    }

    // Peut-etre ne rien faire ?
    if (etat==ETATZERO) {
	return ; 
    }

    assert(etat==ETATQCQ) ; 

    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
	_sxm1[base_r](t[l], base.b[l]) ;
    }
    
}
}
