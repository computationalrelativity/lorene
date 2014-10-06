/*
 * Computation of x*Id
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
 */

/*
 *   Copyright (c) 1999-2001 Jerome Novak
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


char valeur_mult_x_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2008/08/27 08:52:55  jl_cornou
 * Added Jacobi(0,2) polynomials case
 *
 * Revision 1.2  2004/11/23 15:17:19  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.5  1999/11/30  12:45:07  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 1.4  1999/11/24  09:34:16  eric
 * Modif commentaires.
 *
 * Revision 1.3  1999/11/23  16:18:04  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 1.2  1999/11/19  09:31:21  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 1.1  1999/11/16  13:37:24  novak
 * Initial revision
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
void _mult_x_pas_prevu(Tbl *, int &) ;
void _mult_x_r_chebp(Tbl *, int &) ;
void _mult_x_r_chebi(Tbl *, int &) ;
void _mult_x_r_chebpim_p(Tbl *, int &) ;
void _mult_x_r_chebpim_i(Tbl *, int &) ;
void _mult_xm1_cheb(Tbl *, int&) ;
void _mult_x_identite (Tbl *, int &) ;
void _mult_x_r_chebpi_p(Tbl *, int &) ;
void _mult_x_r_chebpi_i(Tbl *, int &) ;
void _mult_x_r_jaco02(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::mult_x() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_mult_x != 0x0) {
	return *p_mult_x ;
    }
    
    // ... si, il faut bosser

    p_mult_x = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_mult_x->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_mult_x->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_mult_x->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->mult_x() ;	// calcul 
    
	p_mult_x->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

   // Termine
    return *p_mult_x ;
}



// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::mult_x() {

// Routines de derivation
static void (*_mult_x[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _mult_x[i] = _mult_x_pas_prevu ;
	}
	// Les routines existantes
	_mult_x[R_CHEB >> TRA_R] = _mult_x_identite ;
	_mult_x[R_CHEBU >> TRA_R] = _mult_xm1_cheb ;
	_mult_x[R_CHEBP >> TRA_R] = _mult_x_r_chebp ;
	_mult_x[R_CHEBI >> TRA_R] = _mult_x_r_chebi ;
	_mult_x[R_CHEBPIM_P >> TRA_R] = _mult_x_r_chebpim_p ;
	_mult_x[R_CHEBPIM_I >> TRA_R] = _mult_x_r_chebpim_i ;
	_mult_x[R_CHEBPI_P >> TRA_R] = _mult_x_r_chebpi_p ;
	_mult_x[R_CHEBPI_I >> TRA_R] = _mult_x_r_chebpi_i ;
	_mult_x[R_JACO02 >> TRA_R] = _mult_x_r_jaco02 ;
    }

    // Debut de la routine 

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    int base_r ;
    for (int l=0 ; l<nzone ; l++) {
	base_r = (base.b[l] & MSQ_R) >> TRA_R ;
	assert(t[l] != 0x0) ;
	_mult_x[base_r](t[l], base.b[l]) ;
    }
}
