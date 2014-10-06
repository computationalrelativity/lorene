/*
 *
 * Multiplication by (x-1) in the external compactified domain
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


char valeur_mult_xm1_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/03/09  16:54:23  eric
 * Traitement du cas etat=ETATZERO
 *
 * Revision 2.2  1999/11/30  12:45:18  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.1  1999/10/18  13:42:33  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.0  1999/04/26  16:16:17  phil
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
void _mult_xm1_identite(Tbl*, int&) ;
void _mult_xm1_cheb(Tbl*, int&) ;


void Valeur::mult_xm1_zec() {

    // Peut-etre ne rien faire ?
    if (etat==ETATZERO) {
	return ; 
    }

    assert(etat==ETATQCQ) ; 

    // Calcul des coef.
    coef() ;
    
    // Division par (x-1) dans la ZEC 
    c_cf->mult_xm1_zec() ;
    set_etat_cf_qcq() ;

    base = c_cf->base ; // On remonte la base de sortie au niveau Valeur
    
}


/*
 * Fonction membre de la classe Mtbl_cf pour la multiplication par (x-1) 
 * dans la zone externe compactifiee applique a this
 *
 */

void Mtbl_cf::mult_xm1_zec()	   
{

// Routines de derivation
static void (*_mult_xm1[MAX_BASE])(Tbl *, int&) ;
static int nap = 0 ;

    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	   _mult_xm1[i] = _mult_xm1_identite ;
	}
	
	// Les routines existantes cas UNSURR
	_mult_xm1[R_CHEBU >> TRA_R] = _mult_xm1_cheb ;
    }

    // Peut-etre ne rien faire ?
    if (etat==ETATZERO) {
	return ; 
    }

    assert(etat==ETATQCQ) ; 

    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
	_mult_xm1[base_r](t[l], base.b[l]) ;
    }
    
}
