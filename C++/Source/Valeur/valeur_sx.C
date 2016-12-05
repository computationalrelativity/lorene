/*
 * Computation of 1/x
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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
 * Revision 1.8  2016/12/05 16:18:21  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2015/03/09 10:32:27  j_novak
 * Inclusion of r-Legendre bases.
 *
 * Revision 1.6  2015/03/05 08:49:33  j_novak
 * Implemented operators with Legendre bases.
 *
 * Revision 1.5  2014/10/13 08:53:51  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:24  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2007/12/14 10:19:35  jl_cornou
 * *** empty log message ***
 *
 * Revision 1.2  2004/11/23 15:17:19  m_forot
 * Added the bases for the cases without any equatorial symmetry
 *  (T_COSSIN_C, T_COSSIN_S, T_LEG, R_CHEBPI_P, R_CHEBPI_I).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  1999/11/30  12:45:46  eric
 * Valeur::base est desormais du type Base_val et non plus Base_val*.
 *
 * Revision 2.4  1999/11/23  16:18:55  eric
 * Reorganisation du calcul dans le cas ETATZERO.
 *
 * Revision 2.3  1999/11/19  09:31:50  eric
 * La valeur de retour est desormais const Valeur &.
 *
 * Revision 2.2  1999/10/28  08:00:45  eric
 * Modif commentaires.
 *
 * Revision 2.1  1999/10/18  13:42:58  eric
 * Suppression de l'argument base dans les routines de derivation des mtbl_cf.
 *
 * Revision 2.0  1999/04/26  14:59:31  phil
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
void _sx_pas_prevu(Tbl *, int &) ;
void _sx_r_chebp(Tbl *, int &) ;
void _sx_r_chebi(Tbl *, int &) ;
void _sx_r_legp(Tbl *, int &) ;
void _sx_r_legi(Tbl *, int &) ;
void _sx_r_chebpim_p(Tbl *, int &) ;
void _sx_r_chebpim_i(Tbl *, int &) ;
void _sxm1_cheb(Tbl *, int&) ;
void _sx_identite (Tbl *, int &) ;
void _sx_r_chebpi_p(Tbl *, int &) ;
void _sx_r_chebpi_i(Tbl *, int &) ;
void _sxpun_r_jaco02(Tbl *, int &) ;
void _sx_r_legp(Tbl *, int &) ;
void _sx_r_legi(Tbl *, int &) ;

// Version membre d'un Valeur
// --------------------------

const Valeur& Valeur::sx() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // Peut-etre rien a faire ?
    if (p_sx != 0x0) {
	return *p_sx ;
    }
    
    // ... si, il faut bosser

    p_sx = new Valeur(mg) ;

    if (etat == ETATZERO) {
	p_sx->set_etat_zero() ; 
    }
    else {
	assert(etat == ETATQCQ) ; 
	p_sx->set_etat_cf_qcq() ;
	Mtbl_cf* cfp = p_sx->c_cf ; // Pointeur sur le Mtbl_cf qui vient d'etre
					// cree par le set_etat_cf_qcq()

	// Initialisation de *cfp : recopie des coef. de la fonction
	if (c_cf == 0x0) {
	    coef() ;
	}
	*cfp = *c_cf ;	
 
	cfp->sx() ;	// calcul 
    
	p_sx->base = cfp->base ; // On remonte la base de sortie au niveau Valeur
    }

    // Termine
    return *p_sx ;
}



// Version membre d'un Mtbl_cf
// ---------------------------

void Mtbl_cf::sx() {

// Routines de derivation
static void (*_sx[MAX_BASE])(Tbl *, int &) ;
static int nap = 0 ;

    // Premier appel
    if (nap==0) {
	nap = 1 ;
	for (int i=0 ; i<MAX_BASE ; i++) {
	    _sx[i] = _sx_pas_prevu ;
	}
	// Les routines existantes
	_sx[R_CHEB >> TRA_R] = _sx_identite ;
	_sx[R_CHEBU >> TRA_R] = _sxm1_cheb ;
	_sx[R_CHEBP >> TRA_R] = _sx_r_chebp ;
	_sx[R_CHEBI >> TRA_R] = _sx_r_chebi ;
	_sx[R_LEG >> TRA_R] = _sx_identite ;
	_sx[R_LEGP >> TRA_R] = _sx_r_legp ;
	_sx[R_LEGI >> TRA_R] = _sx_r_legi ;
	_sx[R_CHEBPIM_P >> TRA_R] = _sx_r_chebpim_p ;
	_sx[R_CHEBPIM_I >> TRA_R] = _sx_r_chebpim_i ;
	_sx[R_CHEBPI_P >> TRA_R] = _sx_r_chebpi_p ;
	_sx[R_CHEBPI_I >> TRA_R] = _sx_r_chebpi_i ;
	_sx[R_JACO02 >> TRA_R] = _sxpun_r_jaco02 ;
	_sx[R_LEG >> TRA_R] = _sx_identite ;
	_sx[R_LEGP >> TRA_R] = _sx_r_legp ;
	_sx[R_LEGI >> TRA_R] = _sx_r_legi ;
    }

    //- Debut de la routine -

    // Protection
    assert(etat == ETATQCQ) ;
    
    // Boucle sur les zones
    for (int l=0 ; l<nzone ; l++) {
	int base_r = (base.b[l] & MSQ_R) >> TRA_R ;
	assert(t[l] != 0x0) ;
	_sx[base_r](t[l], base.b[l]) ;
    }
}
}
